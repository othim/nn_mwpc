
'''
    Oliver Thim 2021
'''
import sys
sys.path.append("/home/x_olith/.local/lib/python3.8/site-packages")

import nn_mwpc
import numpy as np 
import time
import matplotlib.pyplot as plt


def Yam_test(obj):
    print('Testing Yam 1S0')
    
    obj.create_new_potential('my_yam_1','Yamaguchi_1S0')
    obj.print_LECs_in_use('my_yam_1')
    obj.print_params_in_use('my_yam_1')

    beta = 1.0
    lam  = 3.0
    obj.set_params_in_potential('my_yam_1',[beta])
    obj.print_param_values('my_yam_1')

    obj.save_potential_decomposition('my_yam_1')

    obj.set_LECs_in_potential('my_yam_1',[lam])
    obj.print_LEC_values('my_yam_1')
# ------------------------------
# --------- MAIN CODE ----------
print("Constructing object and saving potential")
# Settings
# ------------------------------
potential          = "Yamaguchi_1S0"
Jmax               = 2
cutoff             = 10000.0     # MeV
cut_pow            = 4          # This is the power, n,  in the e^(-p/Lambda)^n regularization
sharp_cutoff       = True      # If true the potential is zer to zero for p>Lambda + 300
precompute_pot     = False      # Precompute and store potential
rel_correction     = False      # If relativistic corrections are implemented
num_grid_points    = 60       # Number of momentum grid points
finite_grid        = False      # If finte momentum grid 
inc_weights_in_pot = True      # Include w and p in potential matrix
cut_on_shell       = True       # Implement the cutoff also on on-shell elements
# -----------------------------

obj = nn_mwpc.nn_mwpc_dwb_interface(potential,Jmax,cutoff,cut_pow,sharp_cutoff,\
        precompute_pot,rel_correction,num_grid_points,finite_grid,\
        inc_weights_in_pot,cut_on_shell)
num_chn = obj.get_chn_len()
print(f'Number of channels: {num_chn}')

Yam_test(obj)
