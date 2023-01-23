'''
    This file tests the M-matrix element
    computation agains nn-on-line for the
    Nijmegen1 potential.

    The Nijmegen1 potential is used in partial waves up to Jmax=9 and then
    the OPE up to Jmax=14. This shows a small error of a few percent in some 
    regions.

    Oliver Thim 2022-12
'''


import sys

import nn_mwpc
import numpy as np 
import time
import matplotlib.pyplot as plt


plt.style.use('../../plot_examples/plt_settings.mplstyle')


def M_matrix(obj,E, ang,S,Mo,Mi):
    
    M_el = obj.compute_M_element(E,ang,S,Mo,Mi)
    
    # Same as in Constants.h (2022-12)
    hbarc = 197.326971941683 
    MeVm2_to_mbarn = (hbarc**2)*10.0
    
    M_el_barn_unit = np.sqrt(MeVm2_to_mbarn)*M_el 

    return M_el_barn_unit

def M_matrix_arr(obj,E, ang_arr,S,Mo,Mi):
    M_arr = np.zeros(ang_arr.shape,dtype=complex)
    for i,a in enumerate(ang_arr):
        m = M_matrix(obj,E,a,S,Mo,Mi)
        M_arr[i] = m
    return M_arr


def read_file_ang_M(filename):

    data = np.loadtxt(filename,skiprows=3)
    return data[:,0],np.array(data[:,1] +1j*data[:,2])

def compare(obj,E,S,Mo,Mi,ax,ax2):
    print(f'Comparing for, S={S}, Mo={Mo}, Mi={Mi}')   
    if (S==0):
        filename = 'np_M_ss_' + str(int(E)) + '_nijm1.txt'
    elif (Mo==0 and Mi==0):
        filename = 'np_M_00_' + str(int(E)) + '_nijm1.txt'
    elif (Mo==1 and Mi==0):
        filename = 'np_M_10_' + str(int(E)) + '_nijm1.txt'
    elif (Mo==0 and Mi==1):
        filename = 'np_M_01_' + str(int(E)) + '_nijm1.txt'
    elif (Mo==1 and Mi==1):
        filename = 'np_M_11_' + str(int(E)) + '_nijm1.txt'
    else:
        print('Error, in quantum numbers')

    print(f'Reading file... + {filename}')
    ang_arr,M_nn = read_file_ang_M(DATA_DIR + filename)
    
    ind = np.arange(0,178)
    ang_arr = np.take(ang_arr,ind)
    M_nn = np.take(M_nn,ind)
    print('Computing M-matrix elements...')
    M_code = M_matrix_arr(obj,E,ang_arr,S,Mo,Mi)
    print('Computing max rel. difference...')
    #print(M_nn)
    #print(M_code)
    max_rel_real = np.max(np.abs(M_code.real-M_nn.real)/np.abs(M_nn.real))
    max_rel_imag = np.max(np.abs(M_code.imag-M_nn.imag)/np.abs(M_nn.imag))
    
    mean_rel_real = np.mean(np.abs(M_code.real-M_nn.real)/np.abs(M_nn.real))
    mean_rel_imag = np.mean(np.abs(M_code.imag-M_nn.imag)/np.abs(M_nn.imag))

    ax2.plot(M_nn.real, color='r',alpha=0.8,label='re nn')
    ax2.plot(M_code.real,color='b',alpha=0.8,label='re code')
    ax2.plot(M_nn.imag, color='r',linestyle='--',alpha=0.8,label='im nn')
    ax2.plot(M_code.imag,color='b',linestyle='--',alpha=0.8,label='im code')

    ax.plot(np.abs(M_nn.real-M_code.real), color='C0',alpha=0.8,label='re abs.err')
    ax.plot(np.abs(M_nn.imag-M_code.imag), color='C1',alpha=0.8,label='im abs.err')
    
    ax.set_title(f'$M^{S}_{{{Mo}{Mi}}}$')
    ax.set_xlabel(r'$\theta_{cm}$ [deg]')
    ax.set_ylabel(r'(mb/sr)$^{1/2}$')
    
    ax2.set_title(f'$M^{S}_{{{Mo}{Mi}}}$')
    ax2.set_xlabel(r'$\theta_{cm}$ [deg]')
    ax2.set_ylabel(r'(mb/sr)$^{1/2}$')
    
    return max_rel_real,max_rel_imag,mean_rel_real,mean_rel_imag
 
# ---------------------------------
# MAIN
# ---------------------------------
print("Constructing object and saving potential")

# Settings
# ------------------------------
potential          = "nijmegen1"
Jmax               = 14
cutoff             = 5000.0     # MeV
cut_pow            = 6          # This is the power, n,  in the e^(-p/Lambda)^n regularization
sharp_cutoff       = False      # If true the potential is zer to zero for p>Lambda + 300
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

# Print info
num_chn = obj.get_chn_len()
print(f'Number of channels: {num_chn}')


DATA_DIR = '../data/'
energies = [10,50,200]

with open('out_M_test.txt', 'w') as f:
    print(f'',file=f)
for E in energies:
    fig,ax = plt.subplots(2,3,figsize=(13,10))
    fig.suptitle(f'E={E} MeV')
    ax = ax.reshape(-1)
    
    fig2,ax2 = plt.subplots(2,3,figsize=(13,10))
    fig2.suptitle(f'E={E} MeV')
    ax2 = ax2.reshape(-1)
    i = 0
    obj.solve_LS_ext_pot(E)
    for S in [0,1]:
        for Mo in range(0,S+1):
            for Mi in range(0,S+1):
                diff = compare(obj,E,S,Mo,Mi,ax[i],ax2[i])
                i = i+1
                with open('out_M_test.txt', 'a') as f:
                    print(f'S= {S}, Mo={Mo}, Mi={Mi}',file=f)
                    print(f'E= {E}, max.rel.diff (real,imag)= {diff[:2]}',file=f)
                    print(f'E= {E}, mean.rel.diff (real,imag)= {diff[2:]}',file=f)
    path = f'fig_M_test_{E}.pdf'
    ax[0].legend()
    fig.tight_layout()
    fig.savefig(path, dpi=500,bbox_inches='tight',pad_inches = 0.1)
    
    path2 = f'fig2_M_test_{E}.pdf'
    ax2[0].legend()
    fig2.tight_layout()
    fig2.savefig(path2, dpi=500,bbox_inches='tight',pad_inches = 0.1)

plt.show()









