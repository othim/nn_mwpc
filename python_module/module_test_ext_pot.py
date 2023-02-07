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


# ------------
# Phase shifts
# ------------
def phase_shifts(obj):
    print('\n \nTesting phase shifts \n')
    # Print the LECs that are in use to know in what order to eneter them in 
    # the calls later.
    #obj.print_LECs_in_use()

    # Just some test values
    C1S0 = -0.112927/100.0
    C3S1 = -0.087340/100.0
    gA2  = 1.29*1.29; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = 1.3e-8;
    C3P2 = 0.1e-8;
    
    LECs = [C1S0,C3P0,C3P2,C3S1,gA2]
    

    obj.solve_LS_ext_pot(10.0)
    for chn_number in [0,0,3]:
        LS_term = obj.get_chn_LS_term(chn_number)
        print(f'LS-Term: {LS_term}')
        start = time.time()
        phases = obj.get_saved_phase_shifts(chn_number)
        end = time.time()
        print(f'Phase shifts in Stapp convention in radians: {phases}')
        print(f'Total time: {1e3*(end-start):0.3f} ms')

def M_matrix(obj):
    
    ang = 140.0 # deg
    E   = 100.0 # MeV

    # Time the funtion call
    start = time.time()
    # First solve the LS equation. This saves phase shifts in the object obj.
    # The saved phase shifts can be accessed with obj.get_saved_phase_shifts(chn_number)
    print('Solving LS eq.')
    obj.solve_LS_ext_pot(E)
    end1 = time.time()
    
    # Call the function that computes an observable at a certain angle. This will
    # be computed with the saves phase shifts from the previous call.
    S  = 0
    Mo = 0
    Mi = 0
    M_el = obj.compute_M_element(E,ang,S,Mo,Mi)
    end = time.time()
    hbarc = 197.326971941683 
    Mevm2_to_mbarn = (hbarc**2)*10.0
    print(f'M_el = {M_el} MeV^-1')    
    print(f'M_el = {M_el*np.sqrt(Mevm2_to_mbarn)} mbarn^1/2')    
    print(f'S={S}, Mo={Mo}, Mi={Mi}') 
    print(f'Total time: {1e3*(end-start):0.3f} ms')
    print(f'Time to solve LS: {1e3*(end1-start):0.3f} ms')
# ------------------------------
# --------- MAIN CODE ----------
print("Constructing object and saving potential")
# Settings
# ------------------------------
potential          = "nijmegen1"
Jmax               = 5
cutoff             = 5000.0     # MeV
cut_pow            = 6          # This is the power, n,  in the e^(-p/Lambda)^n regularization
sharp_cutoff       = False      # If true the potential is set to zero for p>Lambda + 300
precompute_pot     = True       # Precompute and store potential
rel_correction     = False      # If relativistic corrections are implemented
num_grid_points    = 120        # Number of momentum grid points
finite_grid        = False      # If finte momentum grid 
inc_weights_in_pot = False      # Include w and p in potential matrix
cut_on_shell       = True       # Implement the cutoff also on on-shell elements
# -----------------------------
obj = nn_mwpc.nn_mwpc_interface(potential,Jmax,cutoff,cut_pow,sharp_cutoff,\
        precompute_pot,rel_correction,num_grid_points,finite_grid,\
        inc_weights_in_pot,cut_on_shell)
num_chn = obj.get_chn_len()
print(f'Number of channels: {num_chn}')
phase_shifts(obj)
M_matrix(obj)
