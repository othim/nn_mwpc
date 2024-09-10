'''
    This program is an example of how to use the python interface to compute
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
def observables(obs):
    print("\n \nTesting observables \n")
    # Print the LECs that are in use to know in what order to eneter them in 
    # the calls later.
    obj.print_LECs_in_use()

    # Just some test values
    C1S0 = -0.09/100.0
    C3S1 = -0.055/100.0
    gA2  = 1.275*1.275; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = -1e-8;
    C3P2 = -0.2e-8;
        
    LECs = [C1S0,C3P0,C3P2,C3S1,gA2]
    # Test to compute many observables at once
    # -----------------------------------------
    print("Many observables at once")
    angles = np.linspace(1,179,50) # in degrees
    angles_l = angles.tolist()
    energies = [10.0,20.0,100.0] # in MeV
    
    # Time the funtion call
    start = time.time()
    obs_vector = obj.compute_observable_l("I 0000",angles_l,energies,LECs)
    end = time.time()

    #print(obs_vector)
    obj.print_LEC_values() # To confirm that the LECs were set correctly
    print(f'Total time: {1e3*(end-start):0.3f} ms')
    print(f'Time per energy: {1e3*(end-start)/len(energies):0.3f} ms \n')
    # -----------------------------------------
    
    # Test to compute one observable at the time
    # ------------------------------------------
    print("One observable at the time \n")
    ang = 10.0 # deg
    E   = 7.43 # MeV

    # Time the funtion call
    start = time.time()
    # First solve the LS equation. This saves phase shifts in the object obj.
    # The saved phase shifts can be accessed with obj.get_saved_phase_shifts(chn_number)
    obj.solve_LS(E,[C1S0,C3P0,C3P2,C3S1,gA2])
    end1 = time.time()
    print("Solved LS")    
    # Call the function that computes an observable at a certain angle. This will
    # be computed with the saves phase shifts from the previous call.
    obs = obj.compute_observable("SGTL",ang)

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
    #C3S1 = -0.087340/100.0
    C3S1 = -0.13/100.0
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
        phases = np.array(phases)
        print(f'Phase shifts in Stapp convention in deg: {phases*180.0/np.pi}')
        print(f'Total time: {1e3*(end-start):0.3f} ms')

def diagonalization(obj):
    print('\n \nTesting diagonalization \n')


    C1S0 = -0.117/100.0
    #C3S1 = -0.108/100.0
    C3S1 = -0.055/100.0
    gA2  = 1.275*1.275; # Note that gA2 = (gA)^2 and are treated as a LEC.
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

def M_matrix(obj):
    print('\n \nTesting diagonalization \n')
    obj.print_LECs_in_use()

    # Just some test values
    C1S0 = 0
    C3S1 = 0
    gA2  = 1.275*1.275; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = 0
    C3P2 = 0
        
    LECs = [C1S0,C3P0,C3P2,C3S1,gA2]
    
    ang = 70.0 # deg
    E   = 200.0 # MeV

    # Time the funtion call
    start = time.time()
    # First solve the LS equation. This saves phase shifts in the object obj.
    # The saved phase shifts can be accessed with obj.get_saved_phase_shifts(chn_number)
    obj.solve_LS(E,[C1S0,C3P0,C3P2,C3S1,gA2])
    
    end1 = time.time()
    # Call the function that computes an observable at a certain angle. This will
    # be computed with the saves phase shifts from the previous call.
    S  = 0
    Mo = 0
    Mi = 0
    M_el = obj.compute_M_element(E,ang,S,Mo,Mi)
    end = time.time()
    
    # Compute comversion factor from MeV^2 to mbarn
    hbarc = 197.326971941683 
    Mevm2_to_mbarn = (hbarc**2)*10.0

    # Print result
    print(f'M_el = {M_el} MeV^-1 sr^{-1/2}')    
    print(f'M_el = {M_el*Mevm2_to_mbarn} (mbarn/sr)^1/2')    
    print(f'S={S}, Mo={Mo}, Mi={Mi}') 
    print(f'Total time: {1e3*(end-start):0.3f} ms')
    print(f'Time to solve LS: {1e3*(end1-start):0.3f} ms')

# ------------------------------
# --------- MAIN CODE ----------
print("Constructing object and saving potential")
# Settings
# ------------------------------
potential          = "MWPC_LO_J"
Jmax               = 2
cutoff             = 700.0     # MeV
cut_pow            = 4          # This is the power, n,  in the e^(-p/Lambda)^n regularization
sharp_cutoff       = True      # If true the potential is zer to zero for p>Lambda + 300
precompute_pot     = True       # Precompute and store potential
rel_correction     = True      # If relativistic corrections are implemented
num_grid_points    = 60       # Number of momentum grid points
finite_grid        = True      # If finte momentum grid 
inc_weights_in_pot = False      # Include w and p in potential matrix
cut_on_shell       = True       # Implement the cutoff also on on-shell elements
# -----------------------------
obj = nn_mwpc.nn_mwpc_interface(potential,Jmax,cutoff,cut_pow,sharp_cutoff,\
        precompute_pot,rel_correction,num_grid_points,finite_grid,\
        inc_weights_in_pot,cut_on_shell)
#num_chn = obj.get_chn_len()
print(f'Number of channels: {num_chn}')
#phase_shifts(obj)
#observables(obj)
#diagonalization(obj)
#M_matrix(obj)
