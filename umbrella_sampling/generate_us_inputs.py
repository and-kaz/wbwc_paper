#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 23:19:40 2017

@author: kazantsev
"""

### Generate .conf inputs for umbrella sampling ###
### using NAMD and WHAM                 ###
## Assumes the already created .pdb frames of the reaction to sample ##
## Needed: table of initial colvars, pdb frame files, initial (general) .conf
##Pipeline:
#1. read the table of initial colvars into a list;
#2. for each FRAME pdb file :
# put its name into .conf, put its number (i) into .conf output,
# put its number into colvarConfig, put its corresponding 1st and 2nd initial colvars
# into its inp. file;
# copy .conf i folder in run/
###
###

import fileinput
import numpy as np
import os
import re

## The following aproach is copied from StackOverflow
## (https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside)

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
###########################
##### set parameters ######
###########################
sysname0='asite_closed'               ### the name of the system
reacname0='GU_wbwc_neb_full'          ### the name of the reaction
usnum='us1'                           ### needed for multiple runs with different kf
pathcv_lambda=300                     ### MANUALLY enter lambda for pathcv (can select optimal from rmsdmatrix)
num_images=25                         ### number of images for pathcv
us_forceconstant=0.5                  ### raw values (kcal/mol) (not divied by width^2)
us_numsteps=1500000                   ### * integration_step = simulation time   (arbitrary, if .conf assume 12h restarts)
qmdirconst=1250                       ### needed to run multiple QM/MM NAMD runs on the cluster (each requires qm folder to i/o qm files)
###keeping the restraints' constant
zres_const = 100        #  [kcal/mol/A^2], physical force constant for z (deviation from the path)
n1h1_const = 100        #  same, but for n1h1 distance
fing_const = 0.04       #  for opening fingers
## conf file
conffileinit = 'asite_closed_example.conf'
##
vmdfileinit = 'vmd0_hlrn.tcl'
###########################       
#############
# FUNCTIONS #
#############
def open_colvars(sysname, reacname):
    cols=np.transpose(np.loadtxt('colvars/cv_%(sysname)s_%(reacname)s.dat' % vars()))
    return cols
#############
#############
def make_dirs(sysname, reacname, add=''):
    "creates run folders for each PDB frame"
    global uspdb
    uspdb=os.listdir('systems/frames/%(sysname)s/%(reacname)s' % vars())
    uspdb.sort(key=natural_keys)
    runname=sysname + '_' + reacname + add
    if not (os.path.exists('run/%(runname)s' % vars())):
        os.mkdir('run/%(runname)s' % vars())
    ##
    for i in range(len(uspdb)):
        try:
            os.mkdir('run/%(runname)s/%(i)s' % vars())
        except:
            print('clear the folder!!!')
    ##
    return runname
#############
#############
######
def write_confs_hb(runname, sysname, conftempl, colvars, kf = None, mask = 'None', restr='off', restempl='', usnum='us1'):
    "for a predefined number of frames (folders),"
    "use a .conf template to write each conf file"
    "if needed, also writes vmd scripts for restraints."
    " update: also rescales width of the restraitns, so the final kf on them is constant"
    if kf == None:
        kf = us_forceconstant
    #
    if type(mask) == str:
        mask = np.full(len(uspdb), 1)
    ##
    colvars_kf=[]       ## write new colvars table, adding kf
    ## rescaling:
    zres_width, fing_width =  np.sqrt(kf/zres_const), np.sqrt(kf/fing_const)
    for i in range(len(uspdb)):
        if mask[i] == 1: 
            pdbframe=uspdb[i]
            colv1=colvars[0][i]
            colv2=colvars[1][i]
            conf=""
            conf0=fileinput.input("conf/%(conftempl)s" % vars())
            conf1=open('run/%(runname)s/%(i)s/us_%(i)s.conf' % vars(), 'w')
            for line in conf0:
                conf+=(str(line.replace("USSYSTEM", str(sysname)).replace("USFRAMEPDB", str(pdbframe)).replace("USFRAMEOUT", str(i)).replace("USQMDIR",str(qmdirconst+i)).replace("USLAMBDA", str(pathcv_lambda)).replace("USNIMAGES", str(num_images)).replace("USFORCECONST", str(kf)).replace("ZWIDTH", str(zres_width)).replace("FINGWIDTH", str(fing_width)).replace("USNUMSTEPS", str(us_numsteps)).replace("COLVAR1", str(colv1)).replace("COLVAR2", str(colv2))))
            conf1.write(conf)
            conf1.close()
            ##
            if restr=='on':
                res0=fileinput.input("conf/%(restempl)s" % vars())
                res=""
                res1=open('run/%(runname)s/%(i)s/vmdres_%(i)s.tcl' % vars(), 'w')
                for line in res0:
                    res+=(str(line.replace("USSYSTEM", str(sysname)).replace("USFRAMEPDB", str(pdbframe)).replace("USFRAMEOUT", str(i))))
                res1.write(res)
                res1.close()
            ##
            fileinput.close()
            ## write colvars table
            colvars_kf.append(['%(usnum)s_%(i)s' % vars(), colv1, colv2, kf])
    ##
    np.savetxt('colvars/cvkf_%(runname)s_%(kf)s.dat' % vars(), colvars_kf, fmt='%s')
#############
#############
###  RUN  ###
#############
############################
### minima and TS regions with separate Kfs
###########################
###
## 0. specify Kfs for minima and not for minima
kf_0 = 0.5
kf_min = 0.05
## 1. open colvars table
colvars0=open_colvars(sysname0, reacname0)
##
# 2. list the pdb dir, get the number of trajs
runname0=make_dirs(sysname0, reacname0, add='_nebhb')
# 3. specify the mask for the minima
minpos = [[0,10], [20,34]]              ## TO ENTER MANUALLY
mininds = np.array([])
for i in minpos:
    mininds = np.array(np.concatenate((mininds, np.array(range(i[0], i[1])))), dtype = int)
##
minmask = np.array(np.zeros(len(uspdb)), dtype = int)
##
minmask[mininds] = 1            ## mask for MINIMA
############
## 4.0 single Kf, no mask
# write_confs_hb(runname0, sysname0, conffileinit , colvars0, restr='off', restempl=vmdfileinit, usnum='us1')
###
# ### 4.1 write the conf (and vmd) files for TS REGION
# write_confs_hb(runname0, sysname0, conffileinit , colvars0, kf = kf_0, mask = np.abs(minmask - 1), restr='on', restempl=vmdfileinit, usnum='us1')
# ###
# # # ## 4.2 write the conf (and vmd) files for MINIMA REGION
write_confs_hb(runname0, sysname0, conffileinit , colvars0, kf = kf_min, mask = minmask, restr='on', restempl=vmdfileinit, usnum='us2')
# ###








        
