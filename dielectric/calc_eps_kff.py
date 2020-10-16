#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 17:30:50 2019

@author: kazantsev
"""

## Permittivity calculation from dipoles


import fileinput
import time
import numpy as np

import scipy.constants as scon
import os
import ntpath

from pylab import *
import matplotlib


import re



homedir='/home/kazantsev/sci/theor/MM/A_site/epsilon/'
os.chdir(homedir)

############################
### CONSTANTS ##############
############################np.transpose(np.loadtxt('data/water_10k_nvt_dcd50_2ns.dat'))
temperature = 298.0               ## [K]

eps0 = scon.epsilon_0             ## dielectric permittivity of vacuum
kbc = scon.k                      ## Boltzmann constant
npi = scon.pi                     ## pi

## conversion factors
debye_to_Cm = 3.33564e-30          ## [Debye] --> [Couloumb x meter]

############################
### FUNCTIONS ##############
############################
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
######
########
## end of stackoverflow pasta ###

def open_dip_file(filename):
    " loads the dipole traj from the file"
    homedirloc=homedir
    dipole_traj=np.transpose(np.loadtxt('%(homedirloc)s/out/dipole/%(filename)s' % vars()))
    return dipole_traj

def molvol_wat(nmols):
    "only for water: use a molecular volume from literature"
    watmol0 = 2.97 * (1e-10)**3                   ## Angstrom**3  -> m((3))
    return nmols * watmol0
########
########
def cell_vol(x,y,z):
    " the cell volume. assumes input in Angstroem"
    cellvol = x * y * z * (1e-10)**3
    return cellvol
########
########
def parse_volmap(volfile, isovalue, vol0, basefold='out/vol/', unit='angstr', plot_hist='no'):
    ".df file parsing"
    " outputs in m^-3 "
    homedirloc=homedir
    ## tries to open with numpy (faster). if fails (num(data)%3 != 0), then do line by line
    try:
        
        volgrid=np.genfromtxt('%(basefold)s%(volfile)s' % vars(), skip_header=8, skip_footer=1)         ## skipping rows are hardcoded. should be constant for the filetype
        volgrid = volgrid.flatten()
    except ValueError:
        tmp=[]
        for line in open('%(basefold)s%(volfile)s' % vars()):
            if len(line.split())<=3:
                for i in range(len(line.split())):
                    tmp.append(line.split()[i])
        volgrid = np.array(tmp, dtype=float)
    ##
    vfilt=volgrid[volgrid>float(isovalue)]
    ##
    if unit == 'angstr':
        conv_coeff = (1e-10)**3 ## Angstrom**3  -> m((3))
    else:
        conv_coeff = 1
    volume = len(vfilt) * vol0**3 * conv_coeff 
    ##
    if plot_hist=='yes':
        plt.hist(vfilt, bins=200)
    return volume
########
########
def parse_vol_traj(folder, isovalue, vol0, unit='angstr', skip = 1):
    " reads a bunch of .dx files, to get a volume trajectory"
    volfold=os.listdir(folder)
    volfold.sort(key=natural_keys)
    volumes = []
    for i in range(len(volfold))[::skip]:
        volumes.append(parse_volmap('%(folder)s/' % vars() + volfold[i], isovalue, vol0, basefold='', unit=unit, plot_hist='no'))
    ##
    if len(volumes) > 0:
        return volumes
    else:
        raise ValueError('The folder is empty!')
########
######## 
def calc_eps_1(volume, dipole_traj, dippart='all'):
    " calculates per. const."
    " Dielectric response of triplex DNA in ionic solution from simulations, 10.1016/S0006-3495(95)80022-6 "
    ## 0. parse the dipole trajectory (assumes output of 1-col file, Debye units)
    ## select end part of the dipole trajectory for analysis
    if type(dippart)=='str':
        dipole_traj=dipole_traj         ## assume all strings mean "all"
    elif type(dippart) == float:
        endlen=round(len(dipole_traj)*float(dippart))       ## analyze only dippart end of the traj
        dipole_traj=dipole_traj[-endlen:-1]
    elif type(dippart) == list:
        startlen, endlen = round(len(dipole_traj)*float(dippart[0])), round(len(dipole_traj)*float(dippart[1]))     ## specific region
        dipole_traj=dipole_traj[int(startlen):int(endlen)]
    else:
        print('something is wrong with dppart option!')
        return
    dipole_traj_sqr=dipole_traj**2
    dip_mean=np.mean(dipole_traj)
    dip_mean_sqr=np.mean(dipole_traj_sqr)
    ###
    denom = 3 * eps0 * kbc * temperature * volume
    numer = debye_to_Cm**2 * (dip_mean_sqr - dip_mean**2)
    #
    alpha = (numer / denom)
    ### Now two calculations of permittivity (see link above):
    ### 1) 1 + alpha, (for surrounding permittivity = infinity)
    perm1 = 1 + (numer / denom)
    ### 2) for sorrounding permittivity = internal permittivity
    perm2 = (3 * alpha + 1 + (9 * alpha**2 + 6 * alpha + 9)**0.5) / 4
    ##
    return [perm1, perm2]
########
########
def calc_list(namelist, gridspace=1.0, isval0=1.0):
    global epslist
    epslist={}
    for i in namelist:
        sysname=i
        testdipfile = 'asite_dip_%(sysname)s.dat' % vars()
        dip1=open_dip_file(testdipfile)
        wat_voltraj=parse_vol_traj('out/vol/%(sysname)s' % vars(), isoval0, gridspace)
        ##
        epslist[i]=calc_eps_1(max(wat_voltraj), dip1, dippart=0.2)[0]

########
########
def calc_dict_asite_rad_3bp(syslist, distlist, repsdict, gridspace=1.0, isoval=1.0, voltype = 'dist'):
    " hardcoded: systype ->  dist -> replicas -> (so 3-layered dict is hardcoded)"
    " the ladst layer is list: 3 base pairs"
    global epsdict
    global dipdict
    folder='6gsk'
    homedirloc=homedir
    epsdict={}
    dipdict={}
    for i in syslist:
        replist = repsdict[i]
        epsdict[i]={}
        dipdict[i]={}
        for d in distlist:
            epsdict[i][d]={}
            dipdict[i][d]={}
            for n in range(1,4):
                n=str(n)
                epsdict[i][d][n]=[]
                dipdict[i][d][n]=[]
                for r in replist:
                    name0 = '%(i)s_%(r)s_%(d)s_%(n)sbp' % vars()
                    name= 'asite_dip_%(name0)s' % vars()
                    dip1=np.transpose(np.loadtxt('%(homedirloc)sout/dipole/%(folder)s/%(name)s.dat' % vars()))
                    dipdict[i][d][n].append(dip1)
                    #####
                    if voltype == 'traj': 
                        try:
                            wat_voltraj=parse_vol_traj('%(homedirloc)sout/vol/%(folder)s/%(name0)s' % vars(), isoval, gridspace, skip = 6)
                            selvol = np.mean(wat_voltraj)
                            print('std/mean in volume: ', np.std(wat_voltraj) / selvol)
                        except ValueError:
                            epsdict[i][d][n].append(np.nan)
                            continue
                    elif voltype == 'dist':
                        dist_m = d * 1e-10      ## angstrom -> m
                        #dist_m = (d + 5) * 1e-10      ## angstrom -> m       ## ions sphere radius
                        selvol = npi * (4/3) * dist_m**3
                    #diplist.append(dip1)
                    dippart=[0.2,0.7]
                    epsdict[i][d][n].append(calc_eps_1(selvol, dip1, dippart=dippart)[-1])
                    print(name0, ' ', r, ': ', epsdict[i][d][n][-1])
                epsdict[i][d][n]=np.array(epsdict[i][d][n])
########
########
def calc_dict_pol(syslist, distlist, repsdict, voltype = 'traj', gridspace=1.0, isoval=1.0):
    " the ladst layer is replicas"
    global epsdict_pol
    global dipdict_pol
    homedirloc=homedir
    epsdict_pol={}
    dipdict_pol={}
    for i in syslist:
        replist = repsdict[i]
        folder=i
        epsdict_pol[i]={}
        dipdict_pol[i] = {}
        for d in distlist:
            epsdict_pol[i][d]=[]
            dipdict_pol[i][d] = []
            for r in replist:
                name0 = '%(i)s_%(r)s_%(d)s' % vars()
                name= 'asite_dip_%(name0)s' % vars()
                dip1=np.transpose(np.loadtxt('%(homedirloc)sout/dipole/%(folder)s/%(name)s.dat' % vars()))
                dipdict_pol[i][d].append(dip1)
                ######
                if voltype == 'traj': 
                    try:
                        wat_voltraj=parse_vol_traj('%(homedirloc)sout/vol/%(folder)s/%(name0)s' % vars(), isoval, gridspace, skip=6)
                        selvol = np.mean(wat_voltraj)
                        print('std/mean in volume: ', np.std(wat_voltraj) / selvol)
                    except ValueError:
                        epsdict_pol[i][d].append(np.nan)
                        continue
                elif voltype == 'dist':
                    dist_m = d * 1e-10      ## angstrom -> m
                    selvol = npi * (4/3) * dist_m**3
                #diplist.append(dip1)
                #print(name0, ' ', r, ' vol: ', np.mean(wat_voltraj), np.std(wat_voltraj))
                dippart=[0.2,0.7]
                epsdict_pol[i][d].append(calc_eps_1(selvol, dip1, dippart=dippart)[-1])             ## [-1] (perm2) is correct (converges to 80 in water)
                print(name0, ' ', r, ': ', epsdict_pol[i][d][-1])
            epsdict_pol[i][d]=np.array(epsdict_pol[i][d])











