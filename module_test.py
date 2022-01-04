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


import nn_mwpc
import numpy as np 
import time

# -----------
# Observables
# -----------

print("Constructing object")
obj = nn_mwpc.nn_mwpc_interface("MWPC_LO_1",25,450.0,True,True)
#obj = nn_mwpc.nn_mwpc_interface("WPC_LO",25,450.0,True,True)

print("Done")

print(obj.print_LECs_in_use())

C1S0 = -0.112927/100.0
C3S1 = -0.087340/100.0
gA2  = 1.289*1.289;
C3P0 = 0.1;
C3P2 = 0.1;


angles = np.linspace(1,179,50)
angles_l = angles.tolist()
energies = [10.0,50.0,200.0]

# Time the funtion call
start = time.time()
obs_vector = obj.compute_observable("I 0000",angles_l,energies,[C1S0,C3P0,C3P2,C3S1,gA2])
#obs_vector = obj.compute_observable("I 0000",angles_l,energies,[C1S0,C3S1,gA2])
end = time.time()

print(obj.print_LEC_values()) # To confirm they are correct

print(f'Total time: {1e3*(end-start):0.3f} ms,  Per energy: {1e3*(end-start)/len(energies):0.3f} ms')

print(obs_vector)

# ------------
# Phase shifts
# ------------
start = time.time()
phase_shift_vec = obj.compute_phase_shift(73,50.0,[C1S0,C3P0,C3P2,C3S1,gA2]);
end = time.time()
print(f'Time: {1e3*(end-start):0.3f} ms') # - uncoup: 2ms, coup 10ms

print(f'Phase shifts: {phase_shift_vec}')
print('Exit OK')
