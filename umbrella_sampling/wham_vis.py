#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 23:19:40 2017

@author: kazantsev
"""

#############
import fileinput
import subprocess
import numpy as np
import os
import ntpath

#from prody import *
import matplotlib
# matplotlib.rcParams['text.usetex'] = True
from pylab import *

import matplotlib.image as mpimg
import matplotlib.patches as patches
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from scipy.stats import linregress
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
from shutil import copyfile
import math
import scipy.optimize
import pyemma

from matplotlib import cm, tri
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib._png import read_png
from matplotlib.cbook import get_sample_data
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

import re
########## set path to the Grossfield's WHAM code
whampath='/home/kazantsev/sci/app/wham/wham/wham'
wham2dpath='/home/kazantsev/sci/app/wham/wham-2d/wham-2d'
steptlen0 = 0.04 ## length of the traj step in ps
########################
def myround(x, base=2):
    return int(base * round(float(x)/base))
################################################
## The following aproaches are copied from StackOverflow
## (https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside)
def atoi(text):
    return int(text) if text.isdigit() else text
#####
def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
####
## end of stackoverflow copypasta
##################
################## 
def format_trajfiles(folder1, folder2, rep, repto):
    """ changes filenames of trajectories and/or copies them to [folder2]
        also deletes repeats based on the time column"""
    try:
        os.mkdir(folder2)
    except:
        pass
    trajs=os.listdir(folder1)
    trajs.sort(key=natural_keys)
    ##
    for i in range(len(trajs)):
        traj_i='%(folder1)s/' % vars() + trajs[i]
        trajrep=trajs[i]                    ## is copying realy needed?
        ## if multiple replacements, do loop
        if type(rep)==list:
            if len(rep)==len(repto):
                for j in range(len(rep)):
                    trajrep=trajrep.replace(rep[j], repto[j])
            else:
                print('replacement arrays have different length')
                return
        else:
            trajrep=trajrep.replace(rep, repto)
        ##
        traj_j='%(folder2)s/' % vars() + trajrep
        ####
        ## now the filename is done, parse the file to remove repeats
        ## also, REMOVE lines with the wrong number of columns
        ## assumes that the first line (AFTER COL NAMES) is fine (acts as reference)
        ####
        lines_save = set()
        traj_i_read = open(traj_i, 'r')
        traj_j_write = open(traj_j, 'w')
        #
        for line in traj_i_read:
            if line not in lines_save:
                if len(list(lines_save)) > 1:
                    if (len(line) == len(list(lines_save)[1])):
                        traj_j_write.write(line)
                        lines_save.add(line)
                    else:
                        continue
                else:
                    traj_j_write.write(line)
                    lines_save.add(line)
        #
        traj_j_write.close()
        traj_i_read.close()
##################
##################
def get_colvar_tables(file):
    "imports the initial colvar table"
    ## at first, try to read the comments (look for labels)
    try:
        colvars_all=np.loadtxt(file, dtype='str')
    except:
        print('something wrong with the file')
        return
    ##
    try:
        colvars_labels0=np.loadtxt(file, comments=None)[0]
    except:
        colvars_labels0=''
        #print('no labels')
    #
    tlen=len(colvars_all)
    print('colvars loaded. number of expected trajectories: ', tlen, '\n')
    ##
    if colvars_labels0 != '':
        colvars_labels=colvars_labels0.strip(' ').strip('#').split(' ').split(',')                          
        if len(colvars_labels)==len(colvars_all):
            for i in range(len(colvars_labels)):
                print('column ', i, ': ', colvars_labels[i])
    return colvars_all
##################
################## 
def tot_time(Lt, P, sign):
    """calculates sum of elements P[i] if P[i] >= Lt"""
    P=np.array(P, dtype=float) # just to make sure
    O=np.zeros(len(P))
    O[P>=Lt]=Lt
    return sign * np.sum(O) # for minization
##################
################## 
def opt_threschold(lengths):
    """optimizes tot_time(Lt)"""
    opt_Lt = scipy.optimize.minimize(tot_time,np.mean(lengths), method='nelder-mead', args=(lengths, -1), options={'maxiter':1000})
    return opt_Lt['x'][0] # minimize() returns an array
##################
################## 
def traj_sel(colvars, folder1, Lt=None, steplen=steptlen0, plot='yes'):
    """ based on colvars table, and length distribution of the trajs,
        selects those that will be included; also finds max.length.
        returns the list of selected trajs AND the optimal (or selected) threschold
        if Lt (threschold) is not scpecified, find the optimal threschold.
        the optimal threschold maximazes total traj.length"""
    texttraj = os.listdir(folder1)
    texttraj.sort(key=natural_keys)
    texttraj = np.array(texttraj)
    #########
    ## 1. initial check of the trajs
    ## traj IDs must be extracted from the traj filenames. assumes all files are *.colvars.traj
    trajs_id=[]         # all traj from folder, just with IDs
    trajs_names=[]      # those that match the colvars, shortened (for comparison and plotting)
    trajs_fullnames=[]  # those that match the colvars (will return it downstream)
    trajs_len=[]        # lengths of the above
    traj_ext='.colvars.traj'
    for i in range(len(texttraj)):
        trajs_id.append(texttraj[i].replace(ussystem+'_', '').replace(traj_ext, ''))
    trajs_id=np.array(trajs_id)
    ## get all trajs listed in colvars that are present in folder1
    for i in range(len(colvars)):
        if any([n == colvars[i][0] for n in trajs_id]):
            try:
                traj_i='%(folder1)s/' % vars() + texttraj[trajs_id == colvars[i][0]][0]  # tmp var
                tmp_len=len(np.transpose(np.loadtxt(traj_i))[0])
                trajs_names.append(trajs_id[trajs_id == colvars[i][0]][0])  
                trajs_fullnames.append(traj_i) 
                trajs_len.append(tmp_len * steplen)   
            except:
                continue
        else:
            continue
    # now estimate the optimal threschold for trajs (or provide it)
    if Lt == None:
        Lt = opt_threschold(trajs_len)
    else:
        print('manual threschold of %(Lt)s was passed' % vars())
    if plot == 'yes':
        plot_traj_lenhist(trajs_len, trajs_names, Lt, distr='yes')
    #
    traj_len, trajs_fullnames = np.array(trajs_len), np.array(trajs_fullnames)
    traj_filtered=trajs_fullnames[traj_len >= Lt]
    # before Lt was user friendly -- in ps; now convert it back to step numbers
    return traj_filtered, Lt / steplen
################## 
##################
def traj_init(trajs, folderto, trajmask):
    """provided by the list of (filtered) trajectories, copies them to folderto/ as .npy files"""
    for i in trajs:
        traj_load=np.transpose(np.loadtxt(i))[trajmask]
        np.save('%(folderto)s/' % vars() + os.path.splitext(os.path.basename(i))[0] + '.npy', traj_load)
################## 
##################
def timesel_trajfiles(folder1, folder2, start, end, seltype='grow'):
    """ selects part os the trajs and copies them to [folder2].
        trajmask was used in the traj_init.
        if (somehow) the traj is too short,
        it will just report and quit"""
    ## 1. don't overwrite the same files!!
    if folder1==folder2:
        print('select different output folder!')
        return
    try:
        os.mkdir(folder2)
    except:
        pass
    trajs=os.listdir(folder1)
    trajs.sort(key=natural_keys)
    ##
    for i in range(len(trajs)):
        traj_i='%(folder1)s/' % vars() + trajs[i]
        traj_load=np.transpose(np.load(traj_i))
        traj_j='%(folder2)s/' % vars() + trajs[i].split('.')[0] + '.dat'        ## hardcoded file extension. the other functions must adapt
        ## check the length
        ## if the traj is too short, save the end batch of it, the same length as (end - start)
        ## it will only happen for loop traj saving (if end != -1)
        if len(traj_load) >= end:
            tstart, tend = start, end
            np.savetxt(traj_j, traj_load[int(tstart):int(tend)])
        else:
            raise IndexError('trajectory is too short!')
            return
        
##################
### WHAM dealing functions
##################       
def write_wham_inp_gen(wtype, colvars, colmask, trajfolder, whamfolder, trajfilt, cv1_width, cv2_width):
    """ writes input for Grossfield's WHAM;
        it is assumed the trajfolder/ contains already filtered trajs
        colvars file's format assumed: 'traj_identifier [colvar columns] kf_global
        colvar mask must be provided to select colvar columns in the colvar file
        reads the list of filteredOUT trajectories, comments them in the wham input file"""
    global kf1
    global cv2_width_loc
    cv2_width_loc = cv2_width
    ## creat a wham input file
    if wtype=='1d':
        us_wham_inp=open(wham1d_data_blank, 'w')
    elif wtype=='2d':
        us_wham_inp=open(wham2d_data_blank, 'w')
    else:
        raise ValueError('can only do 1d or 2d!')
        return
    if ((wtype=='1d') and len(colmask)!=1) or ((wtype=='2d') and len(colmask)!=2):
        raise ValueError('colvar mask does not correspond to the WHAM type')
        return
    ####  part for getting to relate files to colvar table
    trajs_id=[]
    trajs=os.listdir(trajfolder)
    trajs.sort(key=natural_keys)
    traj_ext='.colvars.npy'
    for i in trajs:
        trajs_id.append(i.replace(ussystem+'_', '').replace(traj_ext, ''))
    trajs_id=np.array(trajs_id)
    ##
    filt=''
    for i in trajs_id:
        cv_i=colvars[np.where(colvars==i)[0][0]]
        ##
        traj_file=whamfolder+'/' + ussystem +'_%s.dat' % cv_i[0]          ## corresponds to timesel_trajfiles()
        ## 3.  set the filter
        if any([n == cv_i[0] for n in trajfilt]):
            filt='#'
        else:
            filt=''
        #####
        if wtype=='1d':
            col1 = float(cv_i[colmask[0]])
            kf1 = (float(cv_i[-1]) / cv1_width**2)
            us_wham_inp.write('%s%s %s %s \n' % (filt, traj_file, col1, kf1))
        else:
            col1, col2 = float(cv_i[colmask[0]]), float(cv_i[colmask[1]])
            kf1, kf2 = ((float(cv_i[-1]) / cv1_width**2)), ((float(cv_i[-1]) / cv2_width**2))       
            us_wham_inp.write('%s%s %s %s %s %s \n' % (filt, traj_file, col1, col2, kf1, kf2))
    #print('successfull')
    us_wham_inp.close()
##################
##################
def iter_filter(colvars, ftype, nsteps):
    """iteratively remove trajectories from wham calculations
        either randomly or in between the present trajectories
        expects colvars directly from get_colvar_tables f-n"""
    colnames=np.transpose(colvars)[0]
    if len(colnames) < nsteps:
        print('number of steps exceeds the number of trajectories')
        return
    indlist=np.array(range(len(colnames)), dtype=int)
    fm=[]                                               ## main list of filters
    if ftype=='random':
        fm.append(['blank'])        ## so the first run would be without filters
        np.random.shuffle(indlist)
        for i in range(nsteps):
            fm.append(fm[i]+[str(colnames[indlist[i]])])
    return fm 
##################
##################                       
def run_wham1d_bash(datafile, xmin, xmax, xnbins, mciters=10):
    """ EXTERNALLY runs Grossfield's WHAM-1D"""
    ## starts bash script with wham
    ## non-periodic is hardcoded !!!
    whampathloc, pmffile0loc = whampath, pmffile1d
    temploc, whamtolloc = temp, whamtol       ## local variables for strings    
    wham_line='%(whampathloc)s Px=0 %(xmin)s %(xmax)s %(xnbins)s %(whamtolloc)s %(temploc)s 0 %(datafile)s %(pmffile0loc)s %(mciters)s 1' % vars()
    wham_line=wham_line.split()
    subprocess.check_output(wham_line)
##################
##################            
def run_wham2d_bash(datafile, xmin, xmax, xbins, ymin, ymax, ybins):
    """ EXTERNALLY runs Grossfield's WHAM-2D"""
    ## non-periodic is hardcoded !!!
    whampathloc, pmffile0loc = wham2dpath, pmffile2d
    temploc, whamtolloc = temp, whamtol       ## local variables for strings
    wham_line='%(whampathloc)s Px=0 %(xmin)s %(xmax)s %(xbins)s Py=0 %(ymin)s %(ymax)s %(ybins)s %(whamtolloc)s %(temploc)s 0 %(datafile)s %(pmffile0loc)s 1' % vars()
    wham_line=wham_line.split()
    subprocess.check_output(wham_line)
##################
################## 
### PyEmma function
##################
def pyemma_preinput(colvars, colmask, cv_widths, trajfolder, trajfilt, middist='yes', exclude=None):
    """
    loads trajs, us centers and Kfs into matrices in a dict
    colvar[colmask] should not contain the time column!
    """
    if len(colmask) != len(cv_widths):
        raise ValueError('number of colums in the mask does not correspond to the number of cv widths')
        return
    cv_widths=np.array(cv_widths)
    us_inpdict={
        "trajs":[],
        "trajlen":[],
        "centers":[],
        "kfs":[]
        }
    ####  part for getting to relate files to colvar table
    trajs_id=[]
    trajs=os.listdir(trajfolder)
    trajs.sort(key=natural_keys)
    traj_ext='.colvars.npy'
    for i in trajs:
        trajs_id.append(i.replace(ussystem+'_', '').replace(traj_ext, ''))
    trajs_id=np.array(trajs_id)
    for i in range(len(trajs_id)):
        cv_i=colvars[np.where(colvars==trajs_id[i])[0][0]]   # a line in colvar file, that corresponds to the traj and contains centers and kfs
        if any([n == cv_i[0] for n in trajfilt]):
            continue    # if this traj was filtered, continue to the next
        else:
            pass
        #####
        us_inpdict['trajs'].append(np.transpose(np.load(trajfolder+'/'+trajs[i]))) # loads the trajectory into the dict
        us_inpdict['trajlen'].append(len(us_inpdict['trajs'][i]))
        us_inpdict['centers'].append(np.array(cv_i[colmask], dtype='float'))
        ## to ease the input, unbiased coords are maked with 0 in cv_width
        kfs = np.array(float(cv_i[-1]) / cv_widths**2)
        kfs[kfs==np.inf] = 0
        us_inpdict['kfs'].append(kfs)
    ###
    # now cut them all to the min distance and exclude the start
    min_length = np.min(us_inpdict['trajlen'])
    # print(min_length)
    for i in range(len(us_inpdict['trajs'])):
        us_inpdict['trajs'][i]=us_inpdict['trajs'][i][exclude:min_length]
    ##
    return us_inpdict

##################
### histo and PMF functions
################## 
def parse_pmf(pmffile, regular='yes', writepmf='no', population='no'):
    """ wham2d writes 9999999.000000 to grid cells with no points;
        so this f-n substitutes it with max(realPMF) + some increment;
        so the PMF can be plotted with regular grid matplotlib f-s"""
    global max_real
    if population=='no':
        popind=3
    elif population=='yes':
        popind=4
    else:
        print('select either yes or no for population column pasing')
        return
    pmffile1=np.transpose(np.loadtxt(pmffile))
    pmf=[]
    max_real=max(pmffile1[2][pmffile1[2] < 100])
    pmfincr=10
    ##
    if regular=='yes':
        pmffile1[2][pmffile1[2] == 9999999.000000] = myround(max_real+pmfincr)
        if writepmf == 'yes':
            np.savetxt('pmfout/pmf_reggrid.dat', np.transpose(pmffile1[0:3]))
        return pmffile1[0:popind]
    elif regular=='no':
        for line in np.transpose(pmffile1):
            if (line[2]!=9999999.000000):
                pmf.append(line)
        return np.transpose(pmf)[0:popind]     
    else:
        print('select grid option')
        return
#################
#################
def calc_wbwc_2d(pmfi, space = 's_hb'):
    """ finds local min in wb and in WC regions
        and their positions (MANUALLY PREDEFINED)
        works for [s, (O6H3-N1H1)]"""
    ## boundaries (wobble should always be first!)
    if space == 's_hb':
        wb_bound=[[0.0, 0.3], [0.3, 1.5]]                 
        wc1_bound=[[0.8, 1.1], [0.0, 1.2]]
        wc2_bound=[[0.8, 1.1], [-1.2, 0]]
        borders = [wb_bound, wc1_bound, wc2_bound]
    elif space == 's_o6h3':
        wb_bound=[[0.0, 0.1], [1.5, 2.2]]
        wc1_bound=[[0.8, 1.1], [1.4, 2.0]]
        borders = [wb_bound, wc1_bound]
    ## energies and coords
    ener, pos = [], []
    for i in range(len(borders)):
        e_i=np.min(pmfi[2][np.where(np.logical_and.reduce((pmfi[0] > borders[i][0][0], pmfi[0] <= borders[i][0][1], pmfi[1] > borders[i][1][0], pmfi[1] <= borders[i][1][1])))])
        ener.append(e_i)
        pos.append(np.transpose(pmfi)[np.where(pmfi[2]==ener[i])][0][0:2])              ## first index is due to np.where
    ## dGf (assumed all other except wb are WC)
    dgf = (np.min(ener[1:]) - ener[0])
    ###
    return [dgf, pos, ener]
        
##################
################## 
def calc_nim_dep(filters, rtype, eqsteps, runsteps, changesteps=None, filtstep=1, filter0='', plotfinpmf='no'):
    """"for a precalculated series of filters, 
        and a series of trajectories lenghts,
        calculate 2d array of wb-WC energy convergence"""
    trajlen=totlen        ## from the traj_init f-n (global variable)
    ##
    if changesteps==None:
        changesteps=runsteps
    ##
    tsteps=np.arange(eqsteps+runsteps, trajlen, changesteps)
    ##
    if rtype=='batches':
        eqchange=np.arange(eqsteps, trajlen-runsteps, changesteps)
    elif rtype=='grow':
        eqchange=np.array([eqsteps for i in range(len(tsteps))])
    ####
    wbwc_matrix=[] ## main matrix of wb-wc chnges
    filt_last=[] ## names of the last elements in the filter (for plotting)
    # now iterate over filters and calculate wb-wc changes over traj. parts
    for i in range(len(filters))[::filtstep]:
        if type(filter0)!=list:
            filter0=[filter0]
        filt_i=filter0 + filters[i]         ## filter0 can contain trajectories that should be turned off constantly
        filt_last.append(filt_i[-1])
        wbwc_change=[] ## array of wb-wc convergence
        ##
        for t in range(len(tsteps)):
            timesel_trajfiles('trajnumpy', wham_folder, eqchange[t], tsteps[t], seltype=rtype)
            write_wham2d_inp(colvars_, [1,2], wham_folder, filt_i, run_wham2d='yes')
            wbwc_e, wb_xy, wc_xy = calc_wbwc(parse_pmf(pmffile0))
            wbwc_change.append(wbwc_e)
            ## plot the final pmf?
            if plotfinpmf=='yes':
                if t==len(tsteps)-1:
                    plot_pmf_cont1(parse_pmf(pmffile0), wbwc='yes', save='yes', saveind='__filt_' + str(i) + '(-' + str(filt_last[-1]) + ')_finpmf')
        ##    
        wbwc_matrix.append(wbwc_change)
    return wbwc_matrix
######################################
############# PLOTTING F-S  ##########
######################################
#################
## For the text positioning
left, width = .1, .84
bottom, height = .1, .84
right = left + width
top = bottom + height
#######################
def plot_traj_lenhist(trajlens, trajnames, traj_lt, distr='yes'):
    """plots histogram of the trajectories lengths;
        steplen [fs] is set to plot lengths in ps;
        will produce 2 plots: individual lengths,
        and total distribution"""
    try:
        os.mkdir('trajplot')
    except:
        pass
    ##
    lengths, traj_lt = np.array(trajlens), traj_lt
    tot_trajlen=np.sum(lengths)
    trcolors=np.array(np.zeros(len(lengths)), dtype='str') ## for bar plots
    trcolors[lengths >= traj_lt] = 'green'
    trcolors[lengths < traj_lt] = 'red'
    incl, excl, tot = len(trcolors[trcolors=='green']), len(trcolors[trcolors=='red']), len(trcolors)
    #print('total traj. length: ', '%.2f' % tot_trajlen, ' ps')
    #
    fig=plt.figure(444,figsize=(4.5,10))
    ax=fig.add_subplot(111)
    ax.barh(list(range(len(lengths))), lengths, color=trcolors, edgecolor=None)
    ax.axvline(traj_lt, ls='dashed', color='darkred', lw=1.6, alpha=0.8)
    #ax.set_ylabel('trajs',fontsize=11, labelpad=1.3)
    ax.set_xlabel('length, ps',fontsize=10, labelpad=1.3)
    plt.yticks(list(range(len(trajnames))), trajnames, fontsize=6)
    ax.set_ylim(-1, len(trajnames)+1)
    ax.xaxis.grid(lw=0.8, ls='dashed', color='gray', alpha=0.7)
    
    plt.subplots_adjust(left=0.27,right=0.99,bottom=0.05,top=0.99,wspace=0.05,hspace=0.05)
    fig.savefig('trajplot/trajlengths.png' % vars(), dpi=300)
    ####
    if distr=='yes':
        fig=plt.figure(455,figsize=(4.7,2.4))
        ax=fig.add_subplot(111)
        ax.hist(lengths, 40, facecolor='gray', edgecolor='dimgray', alpha=0.95, histtype='stepfilled', lw=1.1) ## hardcoded !!!
        ax.axvline(traj_lt, ls='dashed', color='darkred', lw=1.5, alpha=0.95)
        ax.text(0.2,top, 'threschold: %.2f' % traj_lt + '\n excluded %(excl)s / %(tot)s' % vars(),fontsize=9, va='top',ha='left',transform=ax.transAxes)
        ax.set_xlabel('length, ps',fontsize=10, labelpad=1.3)
        plt.subplots_adjust(left=0.12,right=0.99,bottom=0.18,top=0.99,wspace=0.05,hspace=0.05)
        fig.savefig('trajplot/length_dist.png' % vars(), dpi=300)
##################
################## 
def plot_sep_trajs(traj_folder, startpart, ptype='scatt', endpart=-1):
    """"iteratively plots each (part of) trajectory over the total trajectory."
        helps finding outlying trajectories, and how well they overlap"""
    global devgrid
    global xdim
    global ydim
    try:
        os.mkdir('trajplot/hist')
    except:
        pass
    eqsteps=startpart
    endtraj=endpart
    trajall=[]
    traj_files=os.listdir(traj_folder)
    traj_files.sort(key=natural_keys)
    for i in range(len(traj_files)):
        traj_i='%(traj_folder)s/' % vars() + traj_files[i]         ## just a name of the file
        traj_load=np.transpose(np.transpose(np.load(traj_i))[int(startpart):int(endpart)])
        trajall.append(traj_load)
    trajall=np.concatenate(trajall, axis=1)
    ## 
    ## plot the total traj
    fig=plt.figure(99,figsize=(6,5))
    ax=fig.add_subplot(111)
     
    if ptype=='scatt':
        ax.scatter(trajall[1], trajall[2], c='lightgray', s=4, alpha=0.6)  ## it is assumed the the formatted trajs always have three columns: time, col1, col2
    elif ptype=='hist2d':
        ax.hist2d(trajall[1], trajall[2], bins=[cv1_nbins, cv2_nbins], norm=LogNorm(), cmap='cividis', alpha=0.75)
    elif ptype=='deviation':
        ## this is for visualizing deviatio from the path
        if len(trajall) < 4:
            print('trajectories do not have enough columns')        ## 0th column is time, the next three must be s, z, (hb ot o6h3 etc)
            plt.close()
            return
        ##
        xdim, ydim = np.mgrid[np.min(trajall[1]):np.max(trajall[1]):50j, np.min(trajall[3]):np.max(trajall[3]):50j]    # bins number and column indices are hardcoded !!!
        #
        devgrid=griddata(np.transpose([trajall[1], trajall[3]]), trajall[2], (xdim, ydim), method='linear')
        devmin, devmax = np.nanmin(devgrid), np.nanmax(devgrid)
        #print(devmin, devmax)
        devlev=np.arange(devmin, devmax, 0.01)
        refcont = ax.contourf(xdim, ydim, devgrid, levels=devlev, cmap=plt.get_cmap('jet'), alpha=1.0)
        plt.subplots_adjust(left=0.12,right=0.80,bottom=0.15,top=0.90,wspace=0.24,hspace=0.20)
        ### ver cbar
        cax=fig.add_axes([0.85, 0.10, 0.035, 0.55])
        cb=plt.colorbar(refcont, cax=cax, orientation='vertical')
        cb.ax.set_title('z', fontsize=12)
        cb.ax.tick_params(labelsize=8, pad=1.2)
        ##
        for call in refcont.collections: 
                call.remove()
    else:
        print('wrong plotting type')
        return
    #ax.set_xlim(cv1_range[0], cv1_range[-1])
    #ax.set_ylim(cv2_range[0], cv2_range[-1])
    ## now a loop for all trajs
    for i in range(len(traj_files)):
        traj_i='%(traj_folder)s/' % vars() + traj_files[i]         ## just a name of the file
        trajname=traj_files[i].split('.')[0]
        traj_load=np.transpose(np.transpose(np.load(traj_i))[int(startpart):int(endpart)])
        ##
        if ptype =='scatt' or ptype=='hist2d':
            fig.suptitle('%(trajname)s, %(eqsteps)s : %(endtraj)s' % vars(), fontsize=10)
            ptraj_i=ax.scatter(traj_load[1], traj_load[2], c='black', s=4, alpha=1.0)
            fig.savefig('trajplot/hist/traj_%(i)s_%(eqsteps)s__%(endtraj)s' % vars())
            ptraj_i.remove()
        elif ptype =='deviation':
            devgrid_i=griddata(np.transpose([traj_load[1], traj_load[3]]), traj_load[2], (xdim, ydim), method='linear')
            maxdev_i = round(np.nanmax(devgrid_i), 3)
            fig.suptitle('%(trajname)s, %(eqsteps)s : %(endtraj)s\n max Z: %(maxdev_i)s' % vars(), fontsize=10)
            cont_i=ax.contour(xdim, ydim, devgrid_i, levels=devlev, linewidths=0.55, colors='k', alpha=1.0, zorder=100)
            contf_i=ax.contourf(xdim, ydim, devgrid_i, levels=devlev, cmap=plt.get_cmap('jet'), alpha=1.0, zorder=100)
            fig.savefig('trajplot/hist/traj_%(i)s_%(eqsteps)s__%(endtraj)s' % vars())
            for call in cont_i.collections: 
                call.remove()
            for call in contf_i.collections: 
                call.remove()
 
##################
################## 
def plot_pmf_cont1(pmf, space='s_hb', wbwc='no', save='no', saveind='', colmap = 'default'):
    """ 2D contour of the selected pmf. written for [s, hb] profiles.
        if plotting minima pos and dG, expects wbwc to be [wbwc_e, [array of 2d-positions]]"""
    ## det levels
    levels1=np.arange(0, math.ceil(max_real), 1)
    fig=plt.figure(11, figsize=(6,4))
    ax=fig.add_subplot(111)
    ##
    if colmap == 'default':
        colmap = colormap
    cont1=ax.tricontourf(pmf[0], pmf[1], pmf[2], levels=levels1, cmap=colmap)
    ax.tricontour(pmf[0], pmf[1], pmf[2], levels=levels1, cmap='Greys', linewidths=0.15, alpha=0.8)
    ##
    ax.set_xlabel(space.split('_')[0],fontsize=11, labelpad=1.3)
    ax.set_ylabel(space.split('_')[1] + r', $\AA$',fontsize=11, labelpad=1.3)
    ax.xaxis.set_tick_params(labelsize=9)
    ax.yaxis.set_tick_params(labelsize=9)
    ax.set_xlim(-0.05, 1.05) # hardcoded !!!
    ## print energies and wb wc positions
    if wbwc!='no':
        ax.text(0.2,0.15, r'$\Delta G_{wc} =$ ' + str(round(wbwc[0], ndigits=2)) + ' kcal/mol',fontsize=11, va='bottom',ha='left',transform=ax.transAxes)
        #ax.text(0.6, 2.1, s=r'$\Delta$G = ' + str(round(wbwc[0], ndigits=2)) + ' kcal/mol', fontsize=9)
        for j in range(len(wbwc[1])):
            ax.scatter(wbwc[1][j][0], wbwc[1][j][1], s=35, facecolor='black', edgecolor='white',  zorder = 100)
    ## colorbar
    cax=fig.add_axes([0.91, 0.22, 0.027, 0.52])
    cb=plt.colorbar(cont1, cax=cax, orientation='vertical')
    cb.ax.set_title('PMF,'+'\n kcal/mol', fontsize=9, y=1.02, x=0.5)
    cb.ax.tick_params(labelsize=9, pad=0.5)
    ##
    plt.subplots_adjust(left=0.11,right=0.87,bottom=0.13,top=0.98,wspace=0.05,hspace=0.05)
    ##
    if save=='yes':
        plt.savefig('pmfplot/' + ussystem + str(saveind) + '_cont1.png' % vars(), format='png', dpi=300)
        fig.clear()  
##################
################## 
def plot_pmf_cont_final(pmf, space='s_hb', wbwc='no', save='no', saveind='', colmap = 'default', path = 'no'):
    """ also 2D contour, but now additionaly with the minimal E path.
        expects paths as [path_wbwc, path_dpt].
        will plot the path on the countour, and in the inset
    """
    ## det levels
    levels1=np.arange(0, math.ceil(max_real), 1)
    fig=plt.figure(11, figsize=(3.15,1.8))
    ax=fig.add_subplot(111)
    ##
    if colmap == 'default':
        colmap = colormap
    cont1=ax.tricontourf(pmf[0], pmf[1], pmf[2], levels=levels1, cmap=colmap)
    ax.tricontour(pmf[0], pmf[1], pmf[2], levels=levels1, cmap='Greys', linewidths=0.15, alpha=0.8)
    ##
    ax.set_xlabel(space.split('_')[0],fontsize=9, labelpad=1.3)
    ax.set_ylabel(space.split('_')[1] + r', $\AA$',fontsize=9, labelpad=1.3)
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    ax.set_xlim(-0.05, 1.05) # hardcoded !!!
    ## print energies and wb wc positions
    if wbwc!='no':
        #ax.text(0.55,top, r'$\Delta$G = ' + str(round(wbwc[0], ndigits=2)) + ' kcal/mol',fontsize=9, va='top',ha='left',transform=ax.transAxes)
        #ax.text(0.6, 2.1, s=r'$\Delta$G = ' + str(round(wbwc[0], ndigits=2)) + ' kcal/mol', fontsize=9)
        for j in range(len(wbwc[1])):
            ax.scatter(wbwc[1][j][0], wbwc[1][j][1], s=15, facecolor='black', edgecolor='white',  zorder = 105)
    #################
    #### PATH ######
    ### path on the top plot
    ## need to set wb to 0
    wb_norm = path[0][2][0]
    if path != 'no':
        path_tot = np.transpose(np.concatenate((np.transpose(path[0]), np.transpose(path[1]))))
        ax.plot(path[0][0], path[0][1], lw = 1.1, ls = 'dashed', color = 'white', zorder = 100)
        ax.plot(path[1][0], path[1][1], lw = 1.1, ls = 'dashed', color = 'lightgray', zorder = 100)
    ########
    ## inset plot with min FE pathe
        ax2=fig.add_axes([0.33, 0.30, 0.25, 0.20])
        #ax2.set_title(r'$\Delta$G(f), kcal/mol', y=0.92, fontsize=10)
        ax2.tick_params(pad=1.2)
        ax2.yaxis.set_tick_params(labelsize=6, which='both', labelright='on')
        ax2.plot(range(len(path_tot[2])), path_tot[2] - wb_norm, lw=1.2, color='black')
        #ax2.scatter(range(len(wbwc['dgf'])), wbwc['dgf'], s=12, marker='o', facecolor='dimgray', edgecolor='black')
        ax2.yaxis.set_ticks_position('both')
        ax2.yaxis.set_tick_params(labelsize=6)
        ax2.xaxis.set_ticks([])
        ax2.yaxis.set_ticks([0, 10, 20, 30])
        ax2.yaxis.grid(ls='dashed', color='lightgray', lw=0.8, alpha=0.9)
        ax2.set_ylim(-5, np.max(path_tot)-wb_norm+3)
        ax2.set_ylabel(r'$\Delta$G,'+'\n kcal/mol', fontsize=7, labelpad=1.0)
        ax2.set_xlabel('path', fontsize=7, labelpad=1.5)
    ##########
    ##########
    ## colorbar
    cax=fig.add_axes([0.88, 0.22, 0.040, 0.52])
    cb=plt.colorbar(cont1, cax=cax, orientation='vertical')
    cb.ax.set_title('PMF,'+'\n kcal/mol', fontsize=8, y=1.02, x=0.5)
    cb.ax.tick_params(labelsize=8, pad=0.5)
    ##
    plt.subplots_adjust(left=0.18,right=0.82,bottom=0.21,top=0.98,wspace=0.05,hspace=0.05)
    ##
    if save=='yes':
        plt.savefig('pmfplot/' + ussystem + str(saveind) + '_cont1_final.png' % vars(), format='png', dpi=300)
        fig.clear()  
#################
#################
def plot_pmf_3d(pmffile, save='no', saveind='', minima = 'no'):
    """ just a 3D plot of the pmf from the pmffile """
    pmf = parse_pmf(pmffile, regular='no')
    ## det levels
    levels1=np.arange(0, math.ceil(max_real), 0.5)
    ##
    plt.style.use('dark_background')
    figt=plt.figure(787, figsize=(7.2,6))
    axt=figt.add_subplot(111, projection='3d')
    axt.xaxis.pane.fill = False
    axt.yaxis.pane.fill = False
    axt.zaxis.pane.fill = False
    #
    axt.set_xlabel('s',fontsize=11, labelpad=2.3)
    axt.set_ylabel(r'hb, $\AA$',fontsize=11, labelpad=2.3)
    axt.set_zlabel(r'$\Delta$G, kcal/mol',fontsize=11, labelpad=2.3)
    axt.tick_params(pad=2.5)
    axt.grid(False)
    scat3d=axt.scatter(pmf[0], pmf[1], pmf[2], s=16, facecolor = 'yellow', alpha=0.95, edgecolor = 'lightgray', lw=0.2)
    if minima =='yes':
        wbwc_out = calc_wbwc_2d(pmf)
        ## find min
        minpos = np.where(wbwc_out[2] == np.min(wbwc_out[2]))[0]
        for j in range(len(wbwc_out[1])):
            axt.scatter(wbwc_out[1][j][0], wbwc_out[1][j][1], wbwc_out[2][j], marker = 'o', s=55, facecolor='aqua', edgecolor='white', lw=0.5)
            if j != minpos:
                axt.plot([wbwc_out[1][j][0],wbwc_out[1][j][0]], [wbwc_out[1][j][1],wbwc_out[1][j][1]], [-1, wbwc_out[2][j]], lw=1.3, ls='dotted', c='w', alpha=0.7)
                axt.text(wbwc_out[1][j][0], wbwc_out[1][j][1], wbwc_out[2][j]/2, s=round(wbwc_out[2][j],1), fontsize=9, c='w', alpha=0.85)
    axt.set_zlim(-1, math.ceil(np.max(pmf[2])))
    #plt.ginput(1) ## allows selecting points
    plt.subplots_adjust(left=0.1,right=0.90,bottom=0.1,top=0.99,wspace=0.22,hspace=0.25)
    ##
    if save=='yes':
        plt.savefig('pmfplot/' + ussystem + str(saveind) + '_3d.png' % vars(), format='png', dpi=300)
        fig.clear()  

##################
################## 
def plot_us_conv_2d(rtype, eqsteps, runsteps, trajlen, changesteps=None, space = 's_hb', plotconv='yes', plotconts='yes'):
    """" plots PMFs contours for a range of bathes of the trajectories
        batches can be the same lengths (rtype='batches') 
        or with increasing lengths (rtype='grow')
        if changesteps is not specified, batches will not overlap
        finally, plot the dG_wc trajectory"""
    ##
    if changesteps==None:
        changesteps=runsteps
    ##
    tsteps=np.arange(eqsteps+runsteps, trajlen, changesteps)
    ##
    if rtype=='batches':
        eqchange=np.arange(eqsteps, trajlen-runsteps, changesteps)
    elif rtype=='grow':
        eqchange=np.array([eqsteps for i in range(len(tsteps))])
    ##
    wbwc_change=[] ## array to plot convergence
    # get the info for wham
    xmin, xmax = cv1_range
    ymin, ymax = cv2_range
    ##
    for t in range(len(tsteps)):
        timesel_trajfiles('trajnumpy', wham_folder, eqchange[t], tsteps[t], seltype=rtype)
        write_wham_inp_gen('2d', colvars_, [1,2], 'trajnumpy', wham_folder, filtblank, cv1_width=cv1_width0, cv2_width=cv2_width0)
        run_wham2d_bash(wham2d_data_blank, xmin, xmax, cv1_nbins, ymin, ymax, cv2_nbins)
        pmf_tmp = parse_pmf(pmffile2d)
        wbwc_out = calc_wbwc_2d(pmf_tmp, space = space)
        wbwc_change.append(wbwc_out[0])
        if plotconts=='yes':
            plot_pmf_cont1(pmf_tmp, space = space, wbwc=wbwc_out, save='yes', saveind='_'+str(eqchange[t]) + '_' + str(tsteps[t]))     
    ## now plot the convergence
    if plotconv=='yes':
        fig=plt.figure(22, figsize=(4,2.2))
        ax=fig.add_subplot(111)
        plt.plot(tsteps*steptlen0, wbwc_change, lw=1.8, color='black')
        plt.scatter(tsteps*steptlen0, wbwc_change, s=20, facecolor='black')
        #
        plt.xlabel('trajectory length, ps',fontsize=10, labelpad=1.3)
        plt.ylabel(r'$\Delta G_{wc}$,'+'\n kcal/mol',fontsize=10, labelpad=1.3)
        plt.xticks(fontsize=9)
        plt.yticks(fontsize=9)
        plt.grid(color='gray', lw=0.8, alpha=0.7)
        plt.subplots_adjust(left=0.21,right=0.98,bottom=0.17,top=0.98,wspace=0.05,hspace=0.05)
        plt.savefig('pmfplot/' + ussystem + '_conv_' + str(rtype) + '_.png' % vars(), format='png', dpi=300)
        return [tsteps*steptlen0, wbwc_change]
    else:
        return [tsteps*steptlen0, wbwc_change]
    
###########################################################
###########################################################
### set parameters ########
###########################
## system-independent constants
####
Rconst = 1.9858775e-3   # kcal / (mol * K)
temp = 298.0    # K
whamtol = 0.0001
####
## global constants
cv1_width0 = 0.025
cv2_width0 = 0.05
##
cv1_range=np.array([-0.1, 1.1])     # path s
cv2_range=np.array([-1.5, 1.5])     # hb
##
cv1_nbins=int(round(sum(abs(cv1_range)) / cv1_width0))
cv2_nbins=int(round(sum(abs(cv2_range)) / cv2_width0))
##
wham2d_data_blank='wham2d_metadata.text'
##
wham_folder='trajwham'
##
colormap='cividis'
#############################
#### SET THE SYSTEM #########
ussystem='asite_closed'
#############################
pmffile2d='pmf/%(ussystem)s_2d.pmf' % vars()
####
filtblank=[]
###########################
###########################
###    RUN         ########
###########################
plt.style.use('default')
os.chdir('out')
###########
##### 1. prepare the trajectories: put all layers together, remove repeats. may take a while, but only needs once (for final trajs)
format_trajfiles('traj_kf5', 'trajall', ['simulation_id1', 'simulation_name'], ['_us1', ussystem])      ## kf 0.5
format_trajfiles('traj_kf05', 'trajall', ['simulation_id2', 'simulation_name'], ['_us2', ussystem])       ## kf 0.05
#########
###### 2. load the colvar table 
colvars_=get_colvar_tables('cvkf_%(ussystem)s.dat' % vars())
########
####### 3. select trajectories for analysis. exclude the lagging ones. 
#######    the funtions will find the optimal length threschold, if no was specified. the length distributions will be plotted
sel_traj, sel_lt = traj_sel(colvars_,  'trajall')
########
####### 4. convert and save the selected trajs as numpy arrays to speed up calculations of the batches
traj_init(sel_traj, 'trajnumpy', [0,1,3])       ## the mask specifies the columns used for WHAM; in this case, we don't need z column
###################################
###################################
##### 5.  at this point, the trajs are ready for WHAM calculations;
### before that, it is useful to check the historgam of the trajectories
plot_sep_trajs('trajnumpy', 300, ptype='hist2d', endpart=-1)
####
## based on the visualization, you can exclude trajectores by adding them in the filter
# filtblank=['']  ## excluded trajectories
############
###### 6. FINAL WHAM CALCULATIONS ####
### this function will automatically calculate WHAM for specified batches, plot PMF for each batch, and dG_wc trajectory
##############
conv_data = plot_us_conv_2d('batches', 300, 600, int(sel_lt), 600, space = 's_hb', plotconv='yes', plotconts='yes')
############
### 3D pmf 
plot_pmf_3d(pmffile2d, save='no', minima = 'yes')
##########









        
