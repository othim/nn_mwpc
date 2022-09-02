'''
    This program is an example of how to use the python interface to comute
    observables and phase shifts for some given potential models.

    - Imoprt the package and create and object with the potential you want

    - Print the LECs that the potential use, IT IS IN THIS ORDER YOU NEED
      TO ENTER THE VALUES IN THE LIST OF LECS GIVEN TO THE compute_observable
      FUNCTION!

    - To be sure that you have enetered them in the correct order, you can 
      check afterwards that the LEC values are as you expect.

    There are lots of more settings you can change, but that has to be done 
    in the C++ code and then recompile the package. 

    Oliver Thim 2021
'''
import sys
sys.path.append("/home/x_olith/.local/lib/python3.8/site-packages")

import nn_mwpc
import numpy as np 
import time
import matplotlib.pyplot as plt

# -----------
# Observables
# -----------
def compute_observable(obj,T_lab,angle,observable):
    """
        obj   : object to compute from
        T_lab : labe enregy in MeV
        angle : cm-angle in deg
        observable name : Eg. "I 0000" 
    """

    # Use the same test values
    C1S0 = -0.1/100.0
    C3S1 = -0.13/100.0
    gA2  = 1.29*1.29; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = 0;
    C3P2 = 0;
        
    LECs = [C1S0,C3P0,C3P2,C3S1,gA2]
    
    obj.solve_LS(T_lab,[C1S0,C3P0,C3P2,C3S1,gA2])
    obs = obj.compute_observable(observable,angle)

    return obs


# Phase shifts
# ------------
def phase_shifts(obj,T_lab,chn_number):
    """
        Computes phase shifts in the Stapp convention in radians
    """ 
    C1S0 = -0.1/100.0
    C3S1 = -0.13/100.0
    gA2  = 1.29*1.29; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = 0;
    C3P2 = 0;
    
    LECs = [C1S0,C3P0,C3P2,C3S1,gA2]
    
    phases = obj.compute_phase_shift(chn_number,T_lab,LECs)
    return phases

def load_phase_shifts():
    # Read in the check channels
    data_file = 'data/data_gen_corr_Andreas.txt'
    print(f'Reading data: {data_file}')
    Data = np.loadtxt(data_file,skiprows=1)
    
    loc_phase_shifts = {"1S0": 1, "3S1": 2, "E1": 3, "1P1" : 4, "3P1" : 5, \
            "3P0" : 6,"3P2": 7,"E2" : 8, "3D1" : 9, "3F2" : 10}

    return Data, loc_phase_shifts


def run_tests_phase(obj):
    
    # Load phase shifts
    data, loc_phase_shifts = load_phase_shifts()

    # Phase shifts
    T_lab    = np.linspace(1,350,350)
    num_channels = obj.get_chn_len()
    
    phase_shifts_list = []
    for chn_i in range(num_channels):
        phase_shifts_chn = []
        LS_term = obj.get_chn_LS_term(chn_i)
        err = 0
        max_err = np.zeros(4)
        for j,Tl in enumerate(T_lab):
            phases = np.array(phase_shifts(obj,Tl,chn_i))*180/np.pi
            #print(f'{Tl:4.4f}  {phases}') 
            phase_shifts_chn.append(phases)
            
            # If it is an uncoupled channel
            if (len(LS_term)==3):
                err = np.abs(phases[3]-data[j,loc_phase_shifts[LS_term]])/np.abs(phases[3])
                #print(err)
                max_err[3] = np.maximum(max_err[3],err)
            else:
                chns_string = ""
                if (LS_term=="3S-D1"):
                    chns_string = ["3D1","3S1","E1"]
                if (LS_term=="3P-F2"):
                    chns_string = ["3F2","3P2","E2"]

                for k,s in enumerate(chns_string):
                    err = np.abs(phases[k]-data[j,loc_phase_shifts[s]])/np.abs(phases[k])
                    #print(f'{k}  {err}')
                    #print(f'{phases[k]}  {data[j,loc_phase_shifts[s]]}')
                    max_err[k] = np.maximum(max_err[k],err)
        suc = "FAILED"
        if (np.all(max_err)<1e-5):
            suc = "SUCCESS"
        print(f'Channel: {LS_term:<5} | max rel. error (dp, dm,e,d_uncoup): {max_err} | {suc}')


def load_PB():
    data_file = 'data/PB_30_MeV_Andreas_corr.txt'
    print(f'Reading data: {data_file}')
    Data = np.loadtxt(data_file,skiprows=1)
    return Data[:,0], Data[:,3]

def run_tests_PB(obj):
    D_ang, D_PB = load_PB()
    T_lab = 30.0 # MeV
    # Solve LS equation
    # -----------------
    C1S0 = -0.1/100.0
    C3S1 = -0.13/100.0
    gA2  = 1.29*1.29; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = 0;
    C3P2 = 0;
    LECs = [C1S0,C3P0,C3P2,C3S1,gA2]
    obj.solve_LS(T_lab,[C1S0,C3P0,C3P2,C3S1,gA2])
    # ----------------
    err = 0
    max_err = 0
    for i,ang in enumerate(D_ang):
        obs = obj.compute_observable("P n000",ang)
        err = np.abs(obs-D_PB[i])/np.abs(obs)
        max_err = np.maximum(max_err,err)
    suc = "FAILED"
    if (max_err<1e-5):
        suc = "SUCCESS"
    
    print(f'PB max.rel error: {max_err} | {suc}')

# ------------------------------
# --------- MAIN CODE ----------
print("Constructing object and saving potential")

# Same settings as in the test

# Model_name, J_max_chn, cutoff, cut_pow, sharp_cutoff, pre_comp_pot, rel_corr,
# number_of_p_points, finite_grid
obj1 = nn_mwpc.nn_mwpc_interface("MWPC_LO_J",2,500.0,6,False,True,True,120,False)
obj_LO_WPC = nn_mwpc.nn_mwpc_interface("MWPC_LO_1",20,500.0,6,False,True,True,120,False)


# Same settings as I use in the computations

# Model_name, J_max_chn, cutoff, cut_pow, sharp_cutoff, pre_comp_pot, rel_corr,
# number_of_p_points, finite_grid
obj2 = nn_mwpc.nn_mwpc_interface("MWPC_LO_J",2,500.0,6,True,True,True,60,False)

print(f'******************************************************')
print(f'******************************************************')
print(f'********* RUNNING TEST WITH TEST SETTINGS ************')
print(f'******************************************************')
print(f'******************************************************')
print(f'**************** Phase shift test ********************')
run_tests_phase(obj1)
print(f'********************* PB test ************************')
run_tests_PB(obj_LO_WPC)
print(f'******************************************************')
print(f'******************************************************')
print(f'******************************************************')
print(f'****************************************************** \n \n')

'''
print(f'******************************************************')
print(f'******************************************************')
print(f'********* RUNNING TEST WITH RUN SETTINGS *************')
print(f'******************************************************')
print(f'******************************************************')
run_tests_phase(obj2)
print(f'******************************************************')
print(f'******************************************************')
print(f'******************************************************')
print(f'******************************************************')
'''
