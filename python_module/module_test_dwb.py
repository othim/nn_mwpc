
'''
    Oliver Thim 2021
'''
import sys
sys.path.append("/home/x_olith/.local/lib/python3.8/site-packages")

import nn_mwpc
import numpy as np 
import time
import matplotlib.pyplot as plt

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
obj = nn_mwpc.nn_mwpc_dwb_interface(potential,Jmax,cutoff,cut_pow,sharp_cutoff,\
        precompute_pot,rel_correction,num_grid_points,finite_grid,\
        inc_weights_in_pot,cut_on_shell)
num_chn = obj.get_chn_len()
print(f'Number of channels: {num_chn}')
