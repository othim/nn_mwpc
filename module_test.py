import nn_mwpc
import time

obj = nn_mwpc.nn_mwpc_interface("WPC_LO",25,450.0,True,True)


print(obj.print_LECs_in_use())

C1S0 = -0.112927/100.0
C3S1 = -0.087340/100.0


# Time the cuntion call
start = time.time()
obs_vector = obj.compute_observable("I 0000",[1,50,100,150],[50.0],[C1S0,C3S1])
end = time.time()

print(obj.print_LEC_values())

print(f'Total time: {1e6*(end-start)} us')

print(obs_vector)
print('Exit OK')
