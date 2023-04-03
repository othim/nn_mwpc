
'''
    Oliver Thim 2021
'''
import sys
sys.path.append("/home/x_olith/.local/lib/python3.8/site-packages")

import nn_mwpc
import numpy as np 
import time
import matplotlib.pyplot as plt


def add_pot(obj,name,pre_def_name,beta,lam):
    print('\nAddig potential\n-----------------------')
    obj.create_new_potential(name,pre_def_name)
    obj.print_LECs_in_use(name)
    obj.print_params_in_use(name)
    
    # Parameters need to be set prior to saving!
    obj.set_params_in_potential(name,[beta])
    obj.print_param_values(name)

    obj.save_potential_decomposition(name)

    # LECs need to be set after saving!
    obj.set_LECs_in_potential(name,[lam])
    obj.print_LEC_values(name)

    
# ------------------------------
# --------- MAIN CODE ----------
print("Constructing object and saving potential")
# Settings
# ------------------------------
potential          = "Yamaguchi_1S0"
Jmax               = 2
cutoff             = 1000000.0     # MeV
cut_pow            = 10000000          # This is the power, n,  in the e^(-p/Lambda)^n regularization
sharp_cutoff       = False      # If true the potential is zer to zero for p>Lambda + 300
precompute_pot     = False      # Precompute and store potential
rel_correction     = False      # If relativistic corrections are implemented
num_grid_points    = 100       # Number of momentum grid points
finite_grid        = False      # If finte momentum grid 
inc_weights_in_pot = True      # Include w and p in potential matrix
cut_on_shell       = False      # Implement the cutoff also on on-shell elements
# -----------------------------

obj = nn_mwpc.nn_mwpc_dwb_interface(potential,Jmax,cutoff,cut_pow,sharp_cutoff,\
        precompute_pot,rel_correction,num_grid_points,finite_grid,\
        inc_weights_in_pot,cut_on_shell)
num_chn = obj.get_chn_len()
print(f'Number of channels: {num_chn}')

add_pot(obj,'1S0_1','Yamaguchi_1S0',40.0,-10000.0)
add_pot(obj,'1S0_2','Yamaguchi_1S0',40.0,5000.0)


print('\nPotential names\n--------------------')
obj.print_potential_names()

print('\nSolving for the T-matrix\n---------------------')
T = obj.solve_exact_pot_sum_T(1.0,0,'1S0_1','1S0_2')
T = np.array(T)
print(T.shape)
print(f'\n\n{T[-1]}, {(np.abs(T[-1]))**2}\n\n')


for i in range(15):
    T = obj.solve_DWBA_T(1.0,0,i,'1S0_1','1S0_2')
    print(f'i={i}:   {T[-1]}, {(np.abs(T[-1]))**2}')
