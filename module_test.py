import numpy as np # Must be imported first, otherwise segfault occurs....
import nn_mwpc
import time

obj = nn_mwpc.nn_mwpc_interface("WPC_LO",25,450.0,True,True)


print(obj.print_LECs_in_use())

C1S0 = -0.112927/100.0
C3S1 = -0.087340/100.0
gA2  = 1.289*1.289;


angles = np.linspace(1,179,50)
print(angles)
angles_l = angles.tolist()

energies = [10.0,50.0,200.0]
# Time the cuntion call
start = time.time()
#for i in range(100):
obs_vector = obj.compute_observable("I 0000",angles_l,energies,[C1S0,C3S1,gA2])
end = time.time()

print(obj.print_LEC_values()) # To confirm they are correct

print(f'Total time: {1e3*(end-start):0.3f} ms,  Per energy: {1e3*(end-start)/len(energies):0.3f} ms')

print(obs_vector)


print('Exit OK')
