
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
#####     KINETIC MODEL OF    ########
##### CODON-ANTICODON DECODING #######
######################################
### the script is initially designed to setup the model via the graphical interface of Tkinter,
### but later converted to a script to run on a server (that's the reason behind the ugliness of the script);
### you can uncomment the Tkinter part and use the GUI
###########################
#### modules
import time
import numpy as np
import scipy
import scipy.constants as scon
from scipy.integrate import odeint
import tkinter
import tkinter.filedialog
import io
import glob
from pylab import *
import matplotlib
import matplotlib.patches as patches
# from PIL import ImageTk, Image
#from pylab import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backend_bases import key_press_handler
import os
####################
#os.chdir('/home/kazantsev/sci/theor/MM/Ago2/kinetic/')

seed(seed=1)                 ## Set a seed for random numbers
ns=0                         ## Number of solutions (needed for multiple solutions at one instance)
nf=0                         ## Figure number (needed for multiple figures)
############################
### CONSTANTS ##############
############################
temperature = 298.0               ## [K]
Rconst = 1.9858775e-3             # kcal / (mol * K)
RT = temperature * Rconst
eps0 = scon.epsilon_0             ## dielectric permittivity of vacuum
kbc = scon.k                      ## Boltzmann constant
npi = scon.pi                     ## pi
#####################
#####################
#############################
#### DEFAULT PARAMETERS #####
#############################
## addres to the file with the described model
homefold = '/home/kazantsev/kinetic/ribdec'
model0 = '%(homefold)s/sys/ribdec_wbwc.dat' % vars()

R1_0=5e-7                   ## (M) Initiation complexes. Rodnina 2010 (10.1038/emboj.2010.229) (competition assays)
ctRNA0=1.0e-6                 ## (M) cognate tRNA (ternary complexes). Rodnina 2010 (10.1038/emboj.2010.229) (competition assays)
nctRNA0=1.0e-6                ## (M) near-cognate tRNA (ternary complexes). Rodnina 2010 (10.1038/emboj.2010.229) (competition assays) (1 - 3 uM)

degrate = 2e1            ## default degradation rate of tRNAs (doesn't matter. it's just the frequence of 'trnastat')

trange0=4e0               ## Default integration time
tres0=8e6                 ## Default integratio time steps
###
colors=np.tile(['black','red', 'orange', 'green', 'blue', 'violet'], 10)   ### to plot different solution on one Fig.


#############################     
## GLOBAL DICTIONARIES ######
#############################

expdeglist = ['er', 'dr', 'ec', 'dc', 'enc', 'dnc']   ## rate constants for expression and degradation
## as a default, they are all zero (no steady state flow, equilibrium instead)
expdeg = {}
for i in expdeglist:
    expdeg[i] = 0
####
## from lit: all but k2 and q3c/nc are from Rodnina 2017 (http://dx.doi.org/10.1098/rstb.2016.0182)
## k2 and q3c/nc : Ehrenberg 2018 (10.1146/annurev-biophys-060414-034148) (q3nc=3.9e3)
all_kfs = {
    'k1': 1.4e8,
    'q2': 85,
    'k2': 7.2e2,
    'q3c': 2.5e1,
    'q3nc': 3.9e3,
    'k3': 180,
    'q4c': 0.2,
    'q4nc': 140,
    'k4c': 190,
    'k4nc': 0.6,
    'k4wb': 0
    }
####
## dict for wbwc params can be defined either from barries or from rate constants
## at first, initialize all entries (with Nan); then manually enter default entries. then calculate and fill the missing ones
### monomer, open- and closed-center
wbwc_kfs = {
    'mono': {},
    'open': {},
    'closed': {}
    }
for i in wbwc_kfs.keys():
    for j in ['dg', 'ddg', 'hf', 'hr']:
        wbwc_kfs[i][j] = 'nan'
## MANUAL ENTRIES ###
wbwc_kfs['mono']['dg'] = 5
wbwc_kfs['mono']['ddg'] = 10.0
#
wbwc_kfs['open']['dg'] = 10.0
wbwc_kfs['open']['ddg'] = 6.0
#
wbwc_kfs['closed']['dg'] = 3.4
wbwc_kfs['closed']['ddg'] = 6.0
#####
####
## f-s to calculate kfs for in wbwc scheme
def calc_pwc_eq(kf, kr):
    pwc = (kf/kr) / ((kf/kr) + 1)           
    return pwc
##
def calc_kfkr(dg, ddg):
    kf = (scon.k * temperature/scon.h) * np.exp(-ddg/(Rconst*temperature))
    kr = (scon.k * temperature/scon.h) * np.exp(-(ddg - dg)/(Rconst*temperature))
    return kf, kr

def calc_pwc_exp(k4nc, k4c):
    "assumed k4wb = 0"
    pwc = k4nc / k4c
    #pwc = (k4nc / k4c) / ((k4nc/k4c) + (k4c/k4nc))
    return pwc
##
def calc_k4nc(pwc, k4c):
    "reverse to calc_pwc"
    k4nc = pwc * k4c
    return k4nc
##
def calc_qwb(pwc, qncexp, qc):
    "calculate q4nc as a function of pwc"
    " pwcexp is pwc calculate using k4ncexp and k4exp"
    qwb = qc + (qncexp - qc) / (1-pwc)
    return qwb
####
def calc_pwc_t_old(kf, kr, t, wb0 = 1):
    wbeq = (kr/kf) / (kr/kf + 1)
    #wbeq = (kr/kf) 
    pwc = 1 - (wbeq + (wb0 - wbeq) * np.exp(-(kf + kr) * t))
    return pwc
##
def calc_pwc_t(kf, kr, t, wc0 = 1):
    wceq = (kf/kr) / (kf/kr + 1)
    #wbeq = (kr/kf) 
    pwc = wceq + (wc0 - wceq) * np.exp(-(kf + kr) * t)
    return pwc  
###
def calc_kr_eq(kf, pwceq):
    kr = kf/pwceq
    return kr
###
def calc_kr_tau(kf, tau):
    "when 1/(kf+kr) equals the lifetime of the state"
    kr = 1/tau - kf
    return kr
########
def calc_dg_ddg(kfkr):
    dg = - (Rconst*temperature) * np.log(kfkr[0] / kfkr[1])
    #wdict[k]['ddg'] = - (Rconst*temperature) * ((np.log(wdict[k]['hf'])/np.log((scon.k * temperature/scon.h))))
    ddg = - (Rconst*temperature) * np.log(kfkr[0]/(scon.k * temperature/scon.h))
    return [dg, ddg]
    
##############
def vary_rate(ratedict, ratelist, nsteps, minmult = -2, maxmult = 2, space='log'):
    "change rate constants"
    # first, get increments
    if space == 'log':
        incr = np.logspace(minmult, maxmult, nsteps)
    elif space == 'lin':
        incr = np.linspace(minmult, maxmult, nsteps)
    ##
    rdictlist = {}
    for k in ratedict.keys():
        rdictlist[k] = []
    ## loop for each nstep, and change each rate const
    for n in range(nsteps):
        for r in ratedict.keys():
           if any([n == r for n in ratelist]):
               rdictlist[r].append(ratedict[r] * incr[n])
           else:
               rdictlist[r].append(ratedict[r])
    ##
    return rdictlist
######

def fill_wbwc_dict(wdict, fill='auto'):
    """for NaN values in the dict, or for those specified, calculates them from the complementary ones;
    update the full dict of rate constants"""
    global all_kfs
    for k in wdict.keys():
        if ((wdict[k]['dg'] == 'nan') and (wdict[k]['ddg'] == 'nan')) or (fill == 'rates'):
            wdict[k]['dg'] = - (Rconst*temperature) * np.log(wdict[k]['hf']/wdict[k]['hr'])
            #wdict[k]['dg'] = - (Rconst*temperature) * np.log(wdict[k]['hf']/wdict[k]['hr'] / (wdict[k]['hf']/wdict[k]['hr'] + wdict[k]['hr']/wdict[k]['hf'])) #### NORMALIZE!!!!
            #wdict[k]['ddg'] = - (Rconst*temperature) * ((np.log(wdict[k]['hf'])/np.log((scon.k * temperature/scon.h))))
            wdict[k]['ddg'] = - (Rconst*temperature) * np.log(wdict[k]['hf']/(scon.k * temperature/scon.h))
        elif ((wdict[k]['hf'] == 'nan') and (wdict[k]['hr'] == 'nan')) or (fill == 'barriers'):
            wdict[k]['hf'] = (scon.k * temperature/scon.h) * np.exp(-wdict[k]['ddg']/(Rconst*temperature))
            wdict[k]['hr'] = (scon.k * temperature/scon.h) * np.exp(-(wdict[k]['ddg'] - wdict[k]['dg'])/(Rconst*temperature))
    ## add hf and hr as full name rate consts to the full dict
    for s in wdict.keys():
        for k in ['hf', 'hr']:
            all_kfs[k+s[0]] = wdict[s][k]
    return wdict
########################
def set_expdeg():
    global expdeg
    """ checks the variabels for C0, and checkboxes for steady states,
    and fill the dictionary of expression-degradation with corresponding rates"""
    if int(stead_R1) == 1:
        expdeg['dr'] = degrate
        expdeg['er'] = float(Cm0_R1_ent) * degrate
    else:
        expdeg['dr'] = 0
        expdeg['er'] = 0
    if int(stead_Tc) == 1:
        expdeg['dc'] = degrate
        expdeg['ec'] = float(Cm0_ctrna_ent) * degrate
    else:
        expdeg['dc'] = 0
        expdeg['ec'] = 0
    if int(stead_Tnc) == 1:
        expdeg['dnc'] = degrate
        expdeg['enc'] = float(Cm0_nctrna_ent) * degrate
    else:
        expdeg['dnc'] = 0
        expdeg['enc'] = 0

#### 1. Initializing, input parsing and inital conditions ####
    
def init(ret = 'no'):
    # Creates empty arrays and fills X0 arrays.
    # Accepts input from the set_expr
    global X0
    dx = np.array(dxlist)
    X0 = [0 for i in range(len(dx))]
    X0[np.where(dx=='R1')[0][0]] = float(Cm0_R1_ent)
    X0[np.where(dx=='T_c')[0][0]] = float(Cm0_ctrna_ent)
    X0[np.where(dx=='T_nc')[0][0]] = float(Cm0_nctrna_ent)
    if ret == 'yes':
        return X0


def init_vardict():
    ## starts automatically upon running
    # creates an empty dictionary for the ODE output
    global vardict
    global varstatdict
    vardict = {
        'time': []}        
    for i in dxlist:
        vardict[i] = []
    ## manual entries (maybe put as a new func?)
    vardict['rib'] = []         ## sum of all ribosome complexes
    vardict['err'] = []         ## Pw / (Pr + Pw)
    vardict['eta'] = []         ## Pw / Pr
    vardict['pnc'] = []
    vardict['pwc_c4'] = []
    vardict['pnc_wc'] = []
    ###
    ## create the same dict, but for endpoints
    varstatdict = vardict
    ###

#### 2. Preparing for solution, defining the model, integrating and parsing the output ####

def presolve():
    ## parses the solving params.
    ## branches into either single ODE or ODEloop
    ## 1. initiate new (empty) X0, or take if from the last sol
    if int(take_x0) == 1:
        save_eq()
    else:
        init()              
    ###
    set_expdeg()        # set the expression/degradation rates
    t=np.linspace(0, float(trange), int(float(tres)))
    #### propagate time
    if int(take_t0.get()) == 1:
        vardict['time'].append(t + vardict['time'][-1][-1])
    else:
        vardict['time'].append(t)
    ##
    n_ode=(len(dxlist))
    sec=float(trange)
    step=sec/int(float(tres))
    ode(calc_model, X0, t)
################

def presolve_for_rates_set():
    """ opens a new window to select kfs and the range
    for static (endpoints) calculations """
    global typevar
    global entrs
    global Rwind
    Rwind=tkinter.Toplevel()
    Rwind.title('__set_h__')    
    R_name=tkinter.Label(Rwind, height = 2, text='   Enter list of rates to vary', font='arial 9 bold')
    R_name.grid(row=1, column = 1, columnspan = 2)
    ## entry for rates
    klist = tkinter.StringVar(Rwind, value='')
    klistent = tkinter.Entry(Rwind, width = 35, textvariable = klist)
    klistent.grid(row=2, column=1, columnspan = 4)
    ## min/max vals and N steps
    minvallab, maxvallab = tkinter.Label(Rwind, height = 1, text='min *', font='arial 9 bold'), tkinter.Label(Rwind, height = 1, text='max *', font='arial 9 bold')
    minval, maxval = tkinter.StringVar(Rwind, value=-3), tkinter.StringVar(Rwind, value=1)
    minvalent, maxvalent = tkinter.Entry(Rwind, width = 7, textvariable = minval), tkinter.Entry(Rwind, width = 7, textvariable = maxval)
    minvallab.grid(row=3, column=1)
    maxvallab.grid(row=3, column=2, sticky = 'w')
    minvalent.grid(row=4, column=1)
    maxvalent.grid(row=4, column=2, sticky = 'w')
    #
    nsteplab = tkinter.Label(Rwind, height = 1, text='N steps', font='arial 9 bold')
    nsteps = tkinter.StringVar(Rwind, value = 20)
    nstepent = tkinter.Entry(Rwind, width = 6, textvariable = nsteps)
    nsteplab.grid(row=3, column = 3, sticky = 'w')
    nstepent.grid(row=4, column = 3, sticky = 'w')
    ###
    ## space type
    typevar=tkinter.StringVar(Rwind, value='log')
    type_log=tkinter.Radiobutton(Rwind, variable=typevar, value= 'log', text='log', height=2, highlightthickness=0)
    type_lin=tkinter.Radiobutton(Rwind, variable=typevar, value= 'lin', text='lin', height=2, highlightthickness=0)
    type_log.grid(row = 3, column = 4)
    type_lin.grid(row = 4, column = 4)
    ### parse vars and start calculations
    #
    solvbut = tkinter.Button(Rwind, width=5,bg='wheat',text='solve', font='arial 10 bold',command=lambda *args: presove_for_rates_run(all_kfs, klist.get().strip().split(', '), int(nsteps.get()), float(minval.get()), float(maxval.get()), typevar.get()), bd=1)
    solvbut.grid(row=6, column=2, columnspan = 2)



def presove_for_rates_run(ratedict, ratelist, nsteps, minmult, maxmult, space):
    ""
    global all_kfs      ## have to think of better way to do it
    global rdictlist
    ####
    terminal.insert('end', '\n kfs selected: %(ratelist)s \n' % vars(), 'yellow')
    terminal.insert('end', '\n params: %(nsteps)s, %(minmult)s, %(maxmult)s, %(space)s \n' % vars(), 'yellow')
    rdictlist = vary_rate(ratedict, ratelist, nsteps, minmult = minmult, maxmult=maxmult, space=space)
    #####
    for n in range(nsteps):
        ## update the rate consts
        for k in all_kfs.keys():
            all_kfs[k] = rdictlist[k][n]
        if int(take_x0.get()) == 1:
            save_eq()
        else:
            init()              
        ###
        set_expdeg()        # set the expression/degradation rates
        t=np.linspace(0, float(trange.get()), int(float(tres.get())))
        terminal.insert('end', '%(n)s...\n' % vars())
        terminal.update_idletasks()
        ode(calc_model, X0, t, stype='static')
    ####
    terminal.insert('end', 'done \n', 'green')



def odeloop(stype, t):
    ## starts a loop of ODEs for a specific var and a range
    # parsing the variable range and steps
    nsol0=ns+1
    if str(space_type.get())=='lin':
        var_range=np.linspace(float(varmin.get()), float(varmax.get()), int(varres.get()))
    else:
        var_range=np.logspace(float(varmin.get()), float(varmax.get()), int(varres.get()))
    ## branching
    if stype=='apo0':
        for i in var_range:
            X0[0]=i
            ode(risc, X0, t)
    elif stype=='rate':
        for k in var_range:
            for i in range(num_guide):
                for j in range(num_target):
                    if str(k_type.get())=='krmi0':
                        if Fij[i][j]==1:
                            K[i][j][3]=k
                    elif str(k_type.get())=='krmi1':
                        if Fij[i][j]==2:
                            K[i][j][3]=k
                    elif str(k_type.get())=='kmi1':
                        if Fij[i][j]==2:
                            K[i][j][2]=k
                    elif str(k_type.get())=='knc':
                        if Fij[i][j]==4:
                            K[i][j][4]=k
            ode(risc, X0, t)
    # successfull termination     
    nsol=ns
    
###################
####
def read_model(modelfile, dictlist):
    """ opens a file desribing the model, parses it into a dict of d[state]/dt,
    and maps rate constants in strings to the appropriate dicts. Expects a list of strings of dict names
    also, puts the image of the scheme to the tkinter window"""
    global dxdict
    global dxlist
    global import_img
    dxdict, dxlist = {}, []     # the list is needed for fixed ordering
    mod = io.open(modelfile, 'r')
    st = next(mod)
    ### image adress is found
    while 'SCHEME_IMAGE' not in st:
        st = next(mod)
    #image_adress = st.strip().split()[-1]
    #import_img = ImageTk.PhotoImage(Image.open(image_adress).resize((496, 384), Image.ANTIALIAS))
    #scheme.configure(image = import_img)
    ### the file must contain equations for ODE between ***STATES*** and ***END*** statements
    while "***STATES***" not in st:
        st = next(mod)
    #
    while "***END***" not in st:
        st = next(mod)
        try:
            dxdict[st.split('=')[0].strip()] = st.split('=')[1].strip().strip(';')
            dxlist.append(st.split('=')[0].strip())
        except:
            continue
    ## now, add dict names to the equations
    ## also, add state names to the PREDEFINED dict
    for s in dxdict.keys():
        for d in dictlist:
            keys = d + '.keys()'
            for k in eval(keys):
                dxdict[s] = dxdict[s].replace(k, "%(d)s['%(k)s']" % vars())
    ##
    for i in dxdict.keys():
        for j in dxdict.keys():
            if "Xdict['%(j)s']" % vars() not in dxdict[i]:
                dxdict[i] = dxdict[i].replace(j, "Xdict['%(j)s']" % vars())
    modelprint, nstates = os.path.basename(modelfile), len(dxlist)


#########
##
def calc_model(X, t):
    """ this is the function used by odeint. it gets new values propagated by odeint as 1d list,
    parses it into the defined states, then uses these values to calculate new rates by equations given by read_model()"""
    Xdict = {}
    for i in range(len(dxlist)):
        Xdict[dxlist[i]] = X[i]
    dX = []
    for i in range(len(dxlist)):
        dX.append(eval(dxdict[dxlist[i]]))
    return dX
    

def ode (f, x0, t, stype = 'kinetic', dname = 'none'):
    ### The central f-n. Consumes almost all runtime ###
    ## Solves the rate eqs. with the odeint f-n (algorithm is from Fortran)
    global sol0
    sol0=odeint(f, x0, t, rtol=1e-12, atol=1e-12)      ### increased convergence criteria
    sol=np.transpose(sol0)
    if stype == 'kinetic':
        postsolve(sol)
    elif stype == 'static':
        postsolve_stat(sol)
    elif stype == 'manual':
        postsolve_manual(sol, dname)
    elif stype == 'array':
        sol_out = postsolve_array(sol, dname)
        return sol_out

def postsolve(sol):
    ## takes ODEs solution from ode() and parses it into variables.
    ## then appends the vardict with them
    global vardict
    global ns
    ns+=1
    nsol=ns           ### local varible for printing
    sols_list.insert('end','sol %(nsol)s' % vars())
    # parsing the sol into vars
    for i in range(len(np.array(dxlist))):
        vardict[dxlist[i]].append(sol[i])
    ##
    # manual entries
    vardict['rib'].append(np.sum([vardict[i][-1] for i in dxlist if (i[0] != 'P') and (i[0] != 'T')], axis = 0))
    vardict['err'].append(vardict['P_w'][-1] / (vardict['P_r'][-1] + vardict['P_w'][-1]))
    vardict['eta'].append(vardict['P_w'][-1] / vardict['P_r'][-1])
    if 'C4_nc_WC' in vardict.keys():
        vardict['pwc_c4'].append(vardict['C4_nc_WC'][-1] / (vardict['C4_nc_WC'][-1] + vardict['C4_nc_wb'][-1]))
        vardict['pnc_wc'].append(vardict['C4_nc_WC'][-1] / (vardict['C4_nc_WC'][-1] + vardict['C4_c'][-1]))
        vardict['pnc'].append((vardict['C4_nc_WC'][-1] + vardict['C4_nc_wb'][-1]) / (vardict['C4_nc_WC'][-1] + vardict['C4_nc_wb'][-1] + vardict['C4_c'][-1]))
    else:
        vardict['pwc_c4'].append(np.zeros(len(vardict['time'][-1])))
        vardict['pnc_wc'].append(np.zeros(len(vardict['time'][-1])))
        vardict['pnc'].append(vardict['C4_nc'][-1] / (vardict['C4_nc'][-1] + vardict['C4_c'][-1]))
    ## manual entries
    # printing to the terminal
    take_t0_chk.configure(state='normal')
    take_x0_chk.configure(state='normal')
    if str(sol_type.get())=='single':
        terminal.insert('end', '--- ODEs solved. ID = %(nsol)s --- \n' % vars(), 'green')
        plt_button.configure(bg='green', state='normal')
        #save_but.configure(bg='orange', state='normal')
    else:
        terminal.insert('end', '.')
        terminal.update_idletasks()
  
############
def postsolve_stat(sol):
    ## takes ODEs solution from ode() and parses it into variables.
    ## then appends the vardict with them
    global varstatdict

    # parsing the sol into vars
    for i in range(len(np.array(dxlist))):
        varstatdict[dxlist[i]].append(sol[i][-1])
    # manual entries
    varstatdict['rib'].append(np.sum([varstatdict[i][-1] for i in dxlist if (i[0] != 'P') and (i[0] != 'T')], axis = 0))
    varstatdict['err'].append(varstatdict['P_w'][-1]/ (varstatdict['P_r'][-1]+ varstatdict['P_w'][-1]))
    varstatdict['eta'].append(varstatdict['P_w'][-1]/ varstatdict['P_r'][-1])
    if 'C4_nc_WC' in varstatdict.keys():
        varstatdict['pwc_c4'].append(varstatdict['C4_nc_WC'][-1] / (varstatdict['C4_nc_WC'][-1] + varstatdict['C4_nc_wb'][-1]))
        varstatdict['pnc_wc'].append(varstatdict['C4_nc_WC'][-1]/ (varstatdict['C4_nc_WC'][-1] + varstatdict['C4_c'][-1]))
        varstatdict['pnc'].append((varstatdict['C4_nc_WC'][-1] + varstatdict['C4_nc_wb'][-1]) / (varstatdict['C4_nc_WC'][-1] + varstatdict['C4_nc_wb'][-1] + varstatdict['C4_c'][-1]))
    else:
        varstatdict['pwc_c4'].append(0)
        varstatdict['pnc_wc'].append(0)
        varstatdict['pnc'].append(varstatdict['C4_nc'][-1]/ (varstatdict['C4_nc'][-1] + varstatdict['C4_c'][-1]))
    ## manual entries
    # printing to the terminal
    terminal.insert('end', '.')
    terminal.update_idletasks()



def postsolve_manual(sol, dname, statelist = 'auto', printdone = 'yes'):
    """ saves only needed endpoint variables in a specified dictname"""
    global vmdict
    if type(statelist) == str:
        statelist = ['P_r', 'P_w', 'C4_nc_wb', 'C4_nc_WC']

    for i in range(len(np.array(dxlist))):
        if any([n == dxlist[i] for n in statelist]):
            vmdict[dxlist[i]][dname].append(sol[i][-1])
    ##
    vmdict['err'][dname].append(vmdict['P_w'][dname][-1]/ (vmdict['P_r'][dname][-1] + vmdict['P_w'][dname][-1]))
    vmdict['eta'][dname].append(vmdict['P_w'][dname][-1]/ vmdict['P_r'][dname][-1])
    if 'C4_nc_WC' in statelist:
        vmdict['pwc_c4'][dname].append(vmdict['C4_nc_WC'][dname][-1]/ (vmdict['C4_nc_WC'][dname][-1] + vmdict['C4_nc_wb'][dname][-1]))
    ##
    if printdone == 'yes':
        print('finished ' + str(len(vmdict[statelist[0]][dname])) + ' ' + dname)
#####
        
        
def postsolve_array(sol, dname):
    """ saves only needed endpoint variables in a specified dictname;
    assumes dname contains a list of properties to parse/calculate"""
    tmp_sol = {}
    for i in range(len(np.array(dxlist))):
        tmp_sol[dxlist[i]] = sol[i][-1]             ## I anyway need to get all states to properly parse 1d sol
    ##
    sol_out = {}            ## main dict
    for p in dname:
        if any([n == p for n in tmp_sol.keys()]):
            sol_out[p] = tmp_sol[p]                 ## for states, directly take their end-point population
        elif p == 'err':
            sol_out[p] = tmp_sol['P_w'] / (tmp_sol['P_r'] + tmp_sol['P_w'])                      ## error rate
        elif p == 'eta':
            sol_out[p] = tmp_sol['P_w'] / (tmp_sol['P_r'])                      ## error rate
        elif p == 'pwc_c4':
            sol_out[p] = tmp_sol['C4_nc_WC'] / (tmp_sol['C4_nc_WC'] + tmp_sol['C4_nc_wb'])       ## population of WC in C4
        elif p == 'pwc_c3':
            sol_out[p] = tmp_sol['C3_nc_WC'] / (tmp_sol['C3_nc_WC'] + tmp_sol['C3_nc_wb'])       ## population of WC in C3
        elif p == 'pnc_wc':
             sol_out[p] = tmp_sol['C4_nc_WC']/ (tmp_sol['C4_nc_WC'] + tmp_sol['C4_c'])           ## population of C4_WC related to C4_C
        elif p == 'pnc':
            sol_out[p] = (tmp_sol['C4_nc_WC'] + tmp_sol['C4_nc_wb']) / (tmp_sol['C4_nc_WC'] + tmp_sol['C4_nc_wb'] + tmp_sol['C4_c'])
        elif p == 'pnc_c3':
            sol_out[p] = (tmp_sol['C3_nc_WC'] + tmp_sol['C3_nc_wb']) / (tmp_sol['C3_nc_WC'] + tmp_sol['C3_nc_wb'] + tmp_sol['C3_c'])
        else:
            print('no formula for property %(p)s was provided!' % vars())
            sol_out[p] = np.nan
    ##
    return sol_out
#######################
#######################

def manual_calc_kfs(wbwc_list, dec_list, timelen, timesteps):
    """ solves ode without tkinter 
    wbwc_list must be ['statename', ['kf', 'kr']]
    dec_list must be ['rate const', [val1, val2, ..., valN]]
    updates the dicts of the rate constants and solves the model for each supplied value """
    global X0
    global wbwc_kfs
    global all_kfs
    ####
    t = np.linspace(0, timelen, timesteps)
    ##
    ## 1. set new wbwc rates
    wbwc_kfs[wbwc_list[0]]['hf'] = wbwc_list[1][0]
    wbwc_kfs[wbwc_list[0]]['hr'] = wbwc_list[1][1]
    wbwc_kfs = fill_wbwc_dict(wbwc_kfs, fill='rates')
    init() 
    #####
    ## now, loop for each value of a rate const in dec_list
    for k in dec_list[1]:
        all_kfs[dec_list[0]] = k
        kname = dec_list[0] + str(k)
        ##
        #init() 
        set_expdeg()
        ode(calc_model, X0, t, stype='manual', dname = kname)
    ###
    
def manual_calc_array(wbwc_list, dec_kname, dec_list, prop_list, timelen, timesteps, sec_kf = 'none', wbwc_state = 'closed', save = 'yes', savename = 'none'):
    """ solves ode without tkinter 
    wbwc_list must be ['statename', ['kf', 'kr']]
    dec_list must be ['rate const', [val1, val2, ..., valN]]
    updates the dicts of the rate constants and solves the model for each supplied value """
    global wbwc_kfs
    global all_kfs
    ####
    ###########
    ## make folder to write all files
    timestamp = dec_kname + '_' + datetime.datetime.now().strftime("%d_%m_%Y__%H_%M_%S")
    homefold0 = homefold
    ###
    if savename == 'none':
            savename = timestamp
    #####
    os.mkdir('%(homefold0)s/data/%(savename)s' % vars())
    ## create grid for each property
    var_dict = {}
    for p in prop_list:
        var_dict[p] = np.zeros([len(wbwc_list), len(dec_list)])
    ## save the list of wbwc and decoding consts.
    np.savetxt('%(homefold0)s/data/%(savename)s/wbwc_list.txt' % vars(), np.transpose(wbwc_list))
    np.savetxt('%(homefold0)s/data/%(savename)s/dec_list.txt' % vars(), dec_list)
    #####
    t = np.linspace(0, timelen, timesteps)
    X0 = init(ret = 'yes') 
    ##
    for h in range(len(wbwc_list)):
        time_s = time.time()
        ## 1. set new wbwc rates
        wbwc_kfs[wbwc_state]['hf'] = wbwc_list[h][0]
        wbwc_kfs[wbwc_state]['hr'] = wbwc_list[h][1]
        wbwc_kfs = fill_wbwc_dict(wbwc_kfs, fill='rates')
        #####
    ## now, loop for each value of a rate const in dec_list
        for k in range(len(dec_list)):
            all_kfs[dec_kname] = dec_list[k]
            ### modify second rate const. this f-n allows to modify only 2 kfs at a time. if more needed, do it externally;
            ### expects sec_kf = ['rate_const', coeff], where coeff = k2_exp / k1_exp
            if type(sec_kf) != str:
                all_kfs[sec_kf[0]] = dec_list[k] * sec_kf[1]
            ###########
            #init() 
            set_expdeg()
            sol_out = ode(calc_model, X0, t, stype='array', dname = prop_list)
            for p in sol_out.keys():
                var_dict[p][h][k] = sol_out[p]
        time_f = time.time()
        calc_time = round((time_f - time_s)/60, 1)
        print(str(h) + '/' + str(len(wbwc_list)) + ' set finished (%(calc_time)s min.)' % vars())
    ##########
    ## now save the data
    if save == 'yes':
        for p in prop_list:
            np.save('%(homefold0)s/data/%(savename)s/%(p)s.npy' % vars(), var_dict[p])
    return timestamp


def parse_sol_data(folder):
    ## TMP
    homefold0 = homefold
    load_sols = {}
    sol_dir = glob.glob('%(homefold0)s/data/%(folder)s/*.npy' % vars())
    print(sol_dir)
    for i in sol_dir:
        load_sols[os.path.splitext(os.path.basename(i))[0]] = np.load(i)
    return load_sols

def save_eq():
    ## Sets the final values of X from odeint as X0 for the following solution.
    ## COMPATIBLE with set_expr()
    global X0
    X0=sol0[-1]
    terminal.insert('end', '--The final concentrations from the last solution --\n-- were set as the initial concentrations--\n', 'white')

#####################
### PLOTTING ########
#####################

def plot0():
    ## Plots the chosen Y(X) plot for the chosen solutions
    global nf  ## to plot into different windows each time it must remember the order
    plt.style.use('default')
    x=get1
    y=get2
    s=sols_list.curselection()
    nf+=1
    wplot1=tkinter.Toplevel()
    wplot1.title('%(x)s - %(y)s' % vars())
    pl0=figure(nf)
    pl0.suptitle('%(x)s - %(y)s' % vars(), fontsize=12)
    pl1=pl0.add_subplot(111)
    plt.xlabel('%(x)s' % vars(),fontsize=12)
    plt.ylabel('%(y)s'% vars(),fontsize=12)
    space=0
    for i in s:                                  ### If the properties are iterable (for guides and targets) and if not
        if '[' in x and ('[' not in y):
            id1=str(x_id.get()).strip().split(',')
            for a in id1:
                pl1.plot(vardict[x][int(i)][int(a)], vardict[y][int(i)], color=colors[i],linewidth=1.4)
            plt.text(0.1, (0.8-space), ("sol %(i)s" % vars()), color=colors[i], transform=pl1.transAxes)
        
        elif '[' in y and ('[' not in x):
            id2=str(y_id.get()).strip().split(',')
            for b in id2:
                pl1.plot(vardict[x][int(i)], vardict[y][int(i)][int(b)], color=colors[i],linewidth=1.4)
            plt.text(0.1, (0.8-space), ("sol %(i)s" % vars()), color=colors[i], transform=pl1.transAxes)
            
        elif '[' in x and ('[' in y):
            id1=str(x_id.get()).strip().split(',')
            if str(y_id.get()) in Ftypes.keys():
                for a in id1:
                    for b in range(len(vardict['mrna[j]'][int(i)])):
                        if F[int(a)][int(b)]==Ftypes[str(y_id.get())]:
                            pl1.plot(vardict[x][int(i)][int(a)], vardict[y][int(i)][int(b)], color=colors[i],linewidth=1.4)
            else:
                id2=str(y_id.get()).strip().split(',')
                for a in id1:
                    for b in id2:
                        pl1.plot(vardict[x][int(i)][int(a)], vardict[y][int(i)][int(b)], color=colors[i],linewidth=1.4)
            plt.text(0.1, (0.8-space), ("sol %(i)s" % vars()), color=colors[i], transform=pl1.transAxes)
                            
        else:
            pl1.plot(vardict[x][int(i)], vardict[y][int(i)], color=colors[i],linewidth=1.4)
            plt.text(0.1, (0.8-space), ("sol %(i)s" % vars()), color=colors[i], transform=pl1.transAxes)
            
    plt.subplots_adjust(left=0.18,right=0.99,bottom=0.15,top=0.92,wspace=0.05,hspace=0.05)
    space+=(0.5/len(s))
    canvas=FigureCanvasTkAgg(pl0, master=wplot1)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas, wplot1)
    toolbar.update()
    canvas._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

######################################################
##################  TKINTER  #########################
######################################################

##### Create main window ######
# root=tkinter.Tk() 
# #root.tk.call('encoding', 'system', 'unicode')
# #root=tkinter.Toplevel()                    
# root.title('ribdec')
# root.configure(background='black')
# ######
# ##### Title and the scheme
# title=tkinter.Label(root, text='  ', bg='black', fg='lightgray',font='gothic 13 bold')  
# scheme_path='homefold/scheme/scheme_blank.png' % vars()                                                                ### keep this in the same folder as the script
# ## scheme:
# try:
#     img = ImageTk.PhotoImage(Image.open(scheme_path).resize((496, 384), Image.ANTIALIAS))
#     scheme=tkinter.Label(root,image=img)
# except:
#     scheme=tkinter.Label(root, height=16, bg='black', fg='dimgray', text='Please put the file "risc_scheme.png" in the folder with this script\nto display it here. It will help working with the model.', font='arial 9')
#     print('Please put the image "risc_scheme.png" in the folder with this script')

####
### Functions for the Listboxes
def parselist1(event):
    global get1
    ind1=var1_list.curselection()[0]
    get1=var1_list.get(ind1)

def parselist2(event):
    global get2
    ind2=var2_list.curselection()[0]
    get2=var2_list.get(ind2)
    
def parselist3(event):
    global get3
    global ind3
    get3=[]
    ind3=sols_list.curselection()[0]
    get3=sols_list.get(ind3)
###
## F-s for the manual modifications is separate windows
### 1. Expression
def set_k():
    """ Creates a separate toplevel window and shows
    the decoding rate constants, which can be modified;
    does NOT show the wb-wc rates (filtered by h*) """
    global Ewind
    Ewind=tkinter.Toplevel()
    Ewind.title('__set_k__')    
    K_name=tkinter.Label(Ewind, height = 2, text='set rate constants', font='arial 9 bold')
    K_name.grid(row=1, column = 1, columnspan = 2)
    labels = {}
    svars = {}
    entrs = {}
    grid_incr = 2
    for k in all_kfs.keys():
        if k[0]!= 'h':
            labels[k] = tkinter.Label(Ewind, width = 8, text=k)
            svars[k] = tkinter.StringVar(Ewind, value="%.1e" % float(all_kfs[k]))
            entrs[k] = tkinter.Entry(Ewind, width = 10, textvariable = svars[k])
            labels[k].grid(row=grid_incr, column=1, sticky='e')
            entrs[k].grid(row=grid_incr, column=2)
            grid_incr+=1
    ## save the results
    save_k_but=tkinter.Button(Ewind, width=5,bg='wheat',text='save', font='arial 10 bold',command=lambda *args: save_k(svars), bd=1)
    save_k_but.grid(row=grid_incr+1, column=1, columnspan = 2)
####
def set_h():
    global typevar
    global entrs
    """ Creates a separate toplevel window and shows
    the wb-wc reaction rate constants, which can be modified"""
    global Hwind
    Hwind=tkinter.Toplevel()
    Hwind.title('__set_h__')    
    H_name=tkinter.Label(Hwind, height = 2, text='set wb<->WC kinetics', font='arial 9 bold')
    H_name.grid(row=1, column = 1, columnspan = 2)
    ## intry type
    typevar=tkinter.StringVar(Hwind, value='barriers')
    type_barriers=tkinter.Radiobutton(Hwind, variable=typevar, value= 'barriers', command = set_h_types, text=u'\u0394G', height=3, highlightthickness=0)
    type_rates=tkinter.Radiobutton(Hwind, variable=typevar, value= 'rates', command = set_h_types, text='h', height=3, highlightthickness=0)
    type_barriers.grid(row = 1, column = 3)
    type_rates.grid(row = 1, column = 4)
    ### column labels
    col_labels = {}
    col_incr = 2
    for c in [u'\u0394G(wb-wc)', u'\u0394\u0394G(wb-wc)', 'hf', 'hr']:
        col_labels[c] = tkinter.Label(Hwind, width = 12, text=c)
        col_labels[c].grid(row=3, column=col_incr, sticky='e')
        col_incr+= 1
    ######
    labels, svars, entrs = {}, {}, {}
    grid_incr = 4
    for k in wbwc_kfs.keys():
        labels[k] = tkinter.Label(Hwind, width = 8, text=k)
        labels[k].grid(row=grid_incr, column=1, sticky='e')
        svars[k], entrs[k] = {}, {}
        col_incr = 2
        for col in ['dg', 'ddg', 'hf', 'hr']:
            if (col == 'dg' or col =='ddg'):
                svars[k][col] = tkinter.StringVar(Hwind,  value="%.2f" % float(wbwc_kfs[k][col]))
            else:
                svars[k][col] = tkinter.StringVar(Hwind,  value="%.2e" % float(wbwc_kfs[k][col]))
            entrs[k][col] = tkinter.Entry(Hwind, width = 10, textvariable = svars[k][col])
            entrs[k][col].grid(row=grid_incr, column=col_incr)
            col_incr+= 1
        grid_incr+=1
    ## save the results
    save_h_but=tkinter.Button(Hwind, width=5,bg='wheat',text='save', font='arial 10 bold',command=lambda *args: save_h(svars, typevar.get()), bd=1)
    save_h_but.grid(row=grid_incr+1, column=3, columnspan = 2)
    set_h_types()
####
def set_h_types():
    global entrs
    ## enable/disable radiobuttons for a range of input values
    if str(typevar.get())=='rates':
        for k in wbwc_kfs.keys():
            for col in ['hf', 'hr']:
                entrs[k][col].configure(state='normal')
            for col in ['dg', 'ddg']:
                entrs[k][col].configure(state='disabled')

    elif str(typevar.get())=='barriers':
        for k in wbwc_kfs.keys():
            for col in ['hf', 'hr']:
                entrs[k][col].configure(state='disabled')
            for col in ['dg', 'ddg']:
                entrs[k][col].configure(state='normal')
    else:
        for k in wbwc_kfs.keys():
            for col in ['hf', 'hr', 'dg', 'ddg']:
                entrs[k][col].configure(state='disabled')
########
def save_k(svars):
    # saves the changes and calls the init() to parse a new X0
    global all_kfs
    for k in svars.keys():
        all_kfs[k] = float(svars[k].get())
    
    terminal.insert('end', 'rate constants modified\n')
    Ewind.destroy()
##
def save_h(svars, filltype):
    # saves the changes and calls the init() to parse a new X0
    global wbwc_kfs
    for k in svars.keys():
        for c in ['dg', 'ddg', 'hf', 'hr']:
            wbwc_kfs[k][c] = float(svars[k][c].get())
    ## now recalculate the unmodified part
    wbwc_kfs = fill_wbwc_dict(wbwc_kfs, fill=filltype)
    terminal.insert('end', 'wb-wc rate constants modified\n')
    Hwind.destroy()
##########
#### String variables #####
########

### Initial conditions ###
Cm0_R1_ent=R1_0    
Cm0_ctrna_ent=ctRNA0
Cm0_nctrna_ent=nctRNA0
###
# checkboxes for steady state
stead_R1=0
stead_Tc=1
stead_Tnc=1
###
## checkboxes for saving time and X for next soutin
take_t0=0
take_x0=0
### Time ###
trange=trange0
tres=tres0
### Range and resolution of variables for multiple solutions ###

##############
### START ####
##############
## final dict:
wbwc_kfs = fill_wbwc_dict(wbwc_kfs)
#all_kfs = set_wbparams(all_kfs)
all_kfs0 = all_kfs

read_model(model0, ['all_kfs', 'expdeg'])
init_vardict()      # prepare the dict of variables

# tkinter.mainloop()


######################################
######################################
#### MANUAL SOLUTIONS
######################################
######################################
#### 1. RUN
######################################
### set dec range
dec_range = np.logspace(-4, 5, 30)
### properties
prop_list = ['err', 'eta', 'pnc', 'pnc_wc', 'pwc_c4', 'pwc_c3', 'pnc_c3']
### set time
trange, tsteps = 8e3, 4e8
##
k4nc0_coeff = all_kfs['k4nc'] / all_kfs['k4c']
q4nc0_coeff = all_kfs['q4nc'] / all_kfs['q4c']
######################
#### WB-WC PARAMRS
##################
pwc_c3_params_vlow = [7.0, 6] ## barriers
pwc_c3_params_low = [3.4, 6] ## barriers
pwc_c3_params_high = [-2.0, 6] ## barriers
######
## slow
s_bars, sl_bars, sh_bars = [-1, 17.8], [-1, 16.8], [-1, 19.0]           ### slow regimes: [dG_{wc}, dG_{ddag}]
s_rates, sl_rates, sh_rates = [calc_kfkr(s_bars[0], s_bars[1])], [calc_kfkr(sl_bars[0], sl_bars[1])], [calc_kfkr(sh_bars[0], sh_bars[1])]   ## stupid workaround: f-n expects list of lists
####
f_bars = [3.4, 6]                                                   ### fast regime
f_rates = [calc_kfkr(f_bars[0], f_bars[1])]
########
### create dict
wbwc_rdict = {
    'num_s_c3vlow': [s_rates, pwc_c3_params_vlow],
    'num_s_c3low': [s_rates, pwc_c3_params_low],
    'num_sh_c3low': [sh_rates, pwc_c3_params_low],
    'num_sh_c3vlow': [sh_rates, pwc_c3_params_vlow],
    'num_sl_c3vlow': [sl_rates, pwc_c3_params_vlow],
    'num_s_c3high': [s_rates, pwc_c3_params_high],
    'num_f_c3low': [f_rates, pwc_c3_params_low],
    'num_f_c3high': [f_rates, pwc_c3_params_high]
    }
#######
######################
##### RUN ###########
#######################
for r in wbwc_rdict.keys():
    wbwc_kfs['open']['dg'] = wbwc_rdict[r][1][0]
    wbwc_kfs['open']['ddg'] = wbwc_rdict[r][1][1]
    wbwc_kfs = fill_wbwc_dict(wbwc_kfs, fill='barriers')
    manual_calc_array(wbwc_rdict[r][0], 'k4c', dec_range, prop_list, trange, tsteps, sec_kf = 'none', save ='yes', savename = r)
    #manual_calc_array(wbwc_rdict[r][0], 'q4c', dec_range, prop_list, trange, tsteps, sec_kf = ['q4nc', q4nc0_coeff], save ='yes', savename = r)
######################







