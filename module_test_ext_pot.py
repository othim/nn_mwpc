'''
    This program tests the external potentials

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
def observables(obs):
    print("\n \nTesting observables \n")
    print("One observable at the time \n")
    ang = 10.0 # deg
    E   = 200.0 # MeV

    # Time the funtion call
    start = time.time()
    # First solve the LS equation. This saves phase shifts in the object obj.
    # The saved phase shifts can be accessed with obj.get_saved_phase_shifts(chn_number)
    obj.solve_LS_ext_pot(E)
    end1 = time.time()
    print("Solved LS")    
    # Call the function that computes an observable at a certain angle. This will
    # be computed with the saves phase shifts from the previous call.
    obs = obj.compute_observable("A 00kk",ang)
    end = time.time()

    print(f'Observable: {obs}')
    print(f'Total time: {1e3*(end-start):0.3f} ms')
    print(f'Time to solve LS: {1e3*(end1-start):0.3f} ms')

    # ------------------------------------------
#print(obj.print_LEC_values()) # To confirm they are correct


# ------------
# Phase shifts
# ------------
def phase_shifts(obj):
    print('\n \nTesting phase shifts \n')
    # Print the LECs that are in use to know in what order to eneter them in 
    # the calls later.
    obj.print_LECs_in_use()

    # Just some test values
    C1S0 = -0.112927/100.0
    C3S1 = -0.087340/100.0
    gA2  = 1.29*1.29; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = 1.3e-8;
    C3P2 = 0.1e-8;
    
    LECs = [C1S0,C3P0,C3P2,C3S1,gA2]
    
    T_lab = np.linspace(10,50,100) # MeV

    for chn_number in [0,0,3]:
        LS_term = obj.get_chn_LS_term(chn_number)
        print(f'LS-Term: {LS_term}')
        start = time.time()
        phases = obj.compute_phase_shift(chn_number,T_lab[0],LECs)
        end = time.time()
        print(f'Phase shifts in Stapp convention in radians: {phases}')
        print(f'Total time: {1e3*(end-start):0.3f} ms')

def diagonalization(obj):
    print('\n \nTesting diagonalization \n')


    C1S0 = -0.117/100.0
    C3S1 = -0.108/100.0
    gA2  = 1.29*1.29; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = 1.3e-8;
    C3P2 = 0.1e-8;

    LECs = [C1S0,C3P0,C3P2,C3S1,gA2]
    chn_number = 3
    start = time.time()
    eigs = obj.compute_binding_energy(chn_number,LECs)
    end = time.time()
    print(f'Eigenvalues in {obj.get_chn_LS_term(chn_number)} for LECs:')
    obj.print_LEC_values()
    print(np.transpose(eigs))
    print(np.transpose(eigs).shape)
    print(f'Total time: {1e3*(end-start):0.3f} ms')

# ------------------------------
# --------- MAIN CODE ----------
print("Constructing object and saving potential")
obj = nn_mwpc.nn_mwpc_interface("nijmegen1",8,10000.0,6,False,True,True)
num_chn = obj.get_chn_len()
print(f'Number of channels: {num_chn}')
observables(obj)
#diagonalization(obj)
