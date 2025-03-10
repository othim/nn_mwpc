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
def convert_observable_name(obs_name_G):
    '''
        This function converts the observables names
        from the Granada database to the ones used in
        the code.
    '''
    if (obs_name_G == 'DSG'):
        obs_string_C = "I 0000"
    elif(obs_name_G == 'SGT'):
        obs_string_C = "SGT"
    elif(obs_name_G == 'PB'):
        obs_string_C = "P n000"
    elif(obs_name_G == 'CKK'):
        obs_string_C = "A 00mm"
    elif(obs_name_G == 'AYY'):
        obs_string_C = "C nn00"
    elif(obs_name_G == 'AZZ'):
        obs_string_C = "A 00kk"
    elif(obs_name_G == 'DT'):
        obs_string_C = "K 0nn0"
    elif(obs_name_G == 'PT'):
        obs_string_C = "P 0n00"
    elif(obs_name_G == 'SGTL'):
        obs_string_C = "SGTL"
    elif(obs_name_G == 'SGTT'):
        obs_string_C = "SGTT"
    else:
        print("Unknown observable")
        print(obs_name_G)
        obs_string_C = "none"
    return obs_string_C
def observables(obs):
    print("\n \nTesting observables \n")
    # Print the LECs that are in use to know in what order to eneter them in 
    # the calls later.
    obj.print_LECs_in_use()

    # Just some test values
    C1S0 = -0.10768e-2                                                          
    C3S1 = -0.07172e-2  
    gA2  = 1.275*1.275; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = 0#1.3e-8;
    C3P2 = 0#0.1e-8;
        
    LECs = [C1S0,C3S1,C3P0,C3P2]
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

    print(obs_vector)
    #obj.print_LEC_values() # To confirm that the LECs were set correctly
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
    obj.solve_LS(E,LECs)
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

def observables_list(obs):
    print("\n \nTesting observables \n")
    # Print the LECs that are in use to know in what order to eneter them in 
    # the calls later.
    obj.print_LECs_in_use()

    # Just some test values
    C1S0 = -0.10768e-2                                                          
    C3S1 = -0.07172e-2  
    gA2  = 1.275*1.275; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = 0#1.3e-8;
    C3P2 = 0#0.1e-8;
        
    LECs = [C1S0,C3S1,C3P0,C3P2]
    
    obs_arr = ['SGT','DSG','SGTL','SGTT','PB','PT','AYY','AZZ','DT','CKK']
    Tlab_arr = [0.01,3.111264214046822829e+01,1.003411036789297697e+02,2.307715384615384551e+02]
    ang_arr   = [0.0,40,90,130,180]

    print(f'{"obs"} \t {"Tlab (MeV)"} \t {"theta_cm (deg)"} \t {"value (units: mb, mb/sr, 1)"}')
    for Tlab in Tlab_arr:
        # Time the funtion call
        # First solve the LS equation. This saves phase shifts in the object obj.
        # The saved phase shifts can be accessed with obj.get_saved_phase_shifts(chn_number)
        obj.solve_LS(Tlab,LECs)
        for obs in obs_arr:
            for ang in ang_arr:
                value = obj.compute_observable(convert_observable_name(obs),ang)
                print(f'{obs} \t {Tlab} \t {ang} \t {value}')



# ------------
# Phase shifts
# ------------
def phase_shifts(obj):
    print('\n \nTesting phase shifts \n')
    # Print the LECs that are in use to know in what order to eneter them in 
    # the calls later.
    obj.print_LECs_in_use()

    C1S0 = -0.10768e-2                                                          
    C3S1 = -0.07172e-2  
    #gA2  = 1.276*1.276; # Note that gA2 = (gA)^2 and are treated as a LEC.
    gA2  = 0
    C3P0 = 0
    C3P2 = 0
    
    LECs = [C1S0,C3S1,C3P0,C3P2]
    
    #T_lab = np.linspace(10,50,100) # MeV
    #T_lab =1e-2
    T_lab = 1.013311036789297681e+00
    print(f'Tlab={T_lab}')
    for chn_number in [0,1,2,3,4,5,6]:
        LS_term = obj.get_chn_LS_term(chn_number)
        print(f'LS-Term: {LS_term}')
        start = time.time()
        phases = obj.compute_phase_shift(chn_number,T_lab,LECs)
        end = time.time()
        phases = np.array(phases)
        print(f'Phase shifts in Stapp convention in deg: {phases*180.0/np.pi}')
        print(f'Total time: {1e3*(end-start):0.3f} ms')

def diagonalization(obj):
    print('\n \nTesting diagonalization \n')


    C1S0 = -0.10768e-2                                                          
    C3S1 = -0.07172e-2  
    gA2  = 1.275*1.275; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = 0;
    C3P2 = 0;

    LECs = [C1S0,C3S1,C3P0,C3P2]
    chn_number = 3
    start = time.time()
    eigs = obj.compute_binding_energy(chn_number,LECs)
    end = time.time()
    print(f'Eigenvalues in {obj.get_chn_LS_term(chn_number)} for LECs:')
    #obj.print_LEC_values()
    #print(np.transpose(eigs))
    print(f'E = {np.min(eigs)} (MeV)')
    print(np.transpose(eigs).shape)
    print(f'Total time: {1e3*(end-start):0.3f} ms')
    print(f'Carlsson et al. E = −2.211(15) MeV')

def print_all_M(obj):

    Tlab_arr = [0.01,3.111264214046822829e+01,1.003411036789297697e+02,2.307715384615384551e+02]
    ang_arr   = [0.0,40,90,130,180]

    for Tlab in Tlab_arr:
        print(f'\nTlab={Tlab}')
        print(f'{"theta (deg)":15}  {"M_s":31}  {"M_00":31}  {"M_pm":31}  {"M_pp":31}  {"M_p0":31}  {"M_0p":31}')
        for  ang in ang_arr:
            M_s,M_00,M_pm,M_pp,M_p0,M_0p = all_M(obj,Tlab,ang)
            print(f'{ang:<15.8e}  {M_s:<.8e}  {M_00:<.8e}  {M_pm:<.8e}  {M_pp:<.8e}  {M_p0:<.8e}  {M_0p:<.8e}')


def all_M(obj,E,ang):
    if np.abs(ang-90.0)<0.001:
        ang=90.001
    M_pp = M_matrix(obj,E, ang,1,1,1)
    M_00 = M_matrix(obj,E, ang,1,0,0)
    M_pm = M_matrix(obj,E, ang,1,1,-1)
    M_s = M_matrix(obj,E, ang,0,0,0)
    M_p0 = M_matrix(obj,E, ang,1,1,0)
    M_0p = M_matrix(obj,E, ang,1,0,1)

    return M_s,M_00,M_pm,M_pp,M_p0,M_0p

def M_matrix(obj,E,ang,S,Mo,Mi):
    #obj.print_LECs_in_use()

    # Just some test values
    C1S0 = -0.10768e-2                                                          
    C3S1 = -0.07172e-2  
    gA2  = 1.275*1.275; # Note that gA2 = (gA)^2 and are treated as a LEC.
    C3P0 = 0#1.3e-8;
    C3P2 = 0#0.1e-8;
        
    LECs = [C1S0,C3S1,C3P0,C3P2]
    

    # Time the funtion call
    start = time.time()
    # First solve the LS equation. This saves phase shifts in the object obj.
    # The saved phase shifts can be accessed with obj.get_saved_phase_shifts(chn_number)
    obj.solve_LS(E,LECs)
    
    end1 = time.time()
    # Call the function that computes an observable at a certain angle. This will
    # be computed with the saves phase shifts from the previous call.
    M_el = obj.compute_M_element(E,ang,S,Mo,Mi)
    end = time.time()
    
    # Compute comversion factor from MeV^2 to mbarn
    hbarc = 197.326971941683 
    Mevm2_to_mbarn = (hbarc**2)*10.0

    # Print result
    #print(f'M_el = {M_el} MeV^-1 sr^{-1/2}')    
    #print(f'M_el = {M_el*Mevm2_to_mbarn} (mbarn/sr)^1/2')    
    #print(f'S={S}, Mo={Mo}, Mi={Mi}') 
    #print(f'Total time: {1e3*(end-start):0.3f} ms')
    #print(f'Time to solve LS: {1e3*(end1-start):0.3f} ms')

    return M_el

# ------------------------------
# --------- MAIN CODE ----------
print("Constructing object and saving potential")
# Settings
# ------------------------------
potential          = "MWPC_LO_1"
Jmax               = 20
cutoff             = 500.0     # MeV
cut_pow            = 6          # This is the power, n,  in the e^(-p/Lambda)^n regularization
sharp_cutoff       = False     # If true the potential is zer to zero for p>Lambda + 300
precompute_pot     = True       # Precompute and store potential
rel_correction     = True      # If relativistic corrections are implemented
num_grid_points    = 100     # Number of momentum grid points
finite_grid        = True      # If finte momentum grid 
inc_weights_in_pot = False      # Include w and p in potential matrix
cut_on_shell       = True       # Implement the cutoff also on on-shell elements
# -----------------------------

fpi = 92.4
mpi = 138.039
Mp  = 938.2720880259
Mn  = 939.5654203856
inv_fm_to_MeV = 197.3269804

print(f'fpi={fpi}')                                                             
print(f'mpi={mpi}')                                                             
print(f'Mp={Mp}')                                                               
print(f'Mn={Mn}')                                                               
print(f'inv_fm_to_MeV={inv_fm_to_MeV}')                                         
print(f'cutoff={cutoff}')                                                       
print(f'cut_pow={cut_pow}')                                                       
print(f'gA={1.276}')
print(f'Jmax={Jmax}')
print(f'rel_corr={rel_correction}')
print(f'Np={num_grid_points}')

obj = nn_mwpc.nn_mwpc_interface(potential,Jmax,cutoff,cut_pow,sharp_cutoff,\
        precompute_pot,rel_correction,num_grid_points,finite_grid,\
        inc_weights_in_pot,cut_on_shell,fpi,mpi,Mp,Mn,inv_fm_to_MeV)
num_chn = obj.get_chn_len()
print(f'Number of channels: {num_chn}')
phase_shifts(obj)
print_all_M(obj)
#observables(obj)
observables_list(obj)
diagonalization(obj)
#M_matrix(obj)
