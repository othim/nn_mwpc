'''
    This program is an example of how to use the python interface.

    - Imoprt the package and create and object with the potential you want

    - Print the LECs that the potential use, IT IS IN THIS ORDER YOU NEED
      TO ENTER THE VALUES IN THE LIST OF LECS GIVEN TO THE compute_observable
      FUNCTION!

    - To be sure that you have enetered them in the correct order, you can 
      check afterwards that the LEC values are as you expect.

    There are a lote more settings you can change, but that has to be done 
    in the C++ code and recompile the package. 

    Oliver Thim 2021
'''




import nn_mwpc
import numpy as np 
import time


print("Constructing object")
obj = nn_mwpc.nn_mwpc_interface("WPC_LO",25,450.0,True,True)
print("Done")

print(obj.print_LECs_in_use())

C1S0 = -0.112927/100.0
C3S1 = -0.087340/100.0
gA2  = 1.289*1.289;


angles = np.linspace(1,179,50)
angles_l = angles.tolist()
energies = [10.0,50.0,200.0]

# Time the funtion call
start = time.time()
obs_vector = obj.compute_observable("I 0000",angles_l,energies,[C1S0,C3S1,gA2])
end = time.time()

print(obj.print_LEC_values()) # To confirm they are correct

print(f'Total time: {1e3*(end-start):0.3f} ms,  Per energy: {1e3*(end-start)/len(energies):0.3f} ms')

print(obs_vector)

print('Exit OK')
