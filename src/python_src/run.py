#
#
# This file contains the main code for running samplings. First there is some
# loading of data and initialization of parameters. There paramters are just
# related to the sampling and statistical model, the physics model is contained
# in the C++ code that is compiled to the library nn_mwpc.
#
#
# Oliver Thim 2021 -
import sys
import pandas as pd
import numpy as np
import bayes as bayes
import importlib
sys.path.append("/Users/toliver/phd") # Add path to library
import nn_mwpc

# Import experimental data
from nn_database_tools import database_manipulation as dm
fname = '/Users/toliver/phd/nn_database_tools/databases/granada2013_absolute_errors_granada_norms.dat'


# -------------------------------------
# Initialize constants
# -------------------------------------

settings['lambda_cut']  = 450.0
settings['lambda_chi']  = 600.0
settings['c_bar']       = 5.0
settings['width_prior'] = 5.0
settings['J_max_saved'] = 25
settings['rel_corr']    = True

settings['s_ndim']      = 4
settings['s_nwalkers']  = 50
settings['s_nsteps']    = 1000
settings['s_initial']   = [1,1,1,1]

# Construct the model to acess the computation of observables
# The hard coded 'True' is if the potential should be pre-computed
# and stored. This does not affect the program in any way except it running
# faster.
model = nn_mwpc.nn_mwpc_interface("WPC_LO",settings['J_max_saved'], \
        settings['lambda_cut'], True, settings['rel_corr'])

# Print so the LECs are in the correct order
model.print_LECs_in_use();

# Unit conversion from the 'GeV' units to MeV units that are accepted in 
# the C++ code. It is very important that the order here corresponds to 
# the order of the LECs as printed by the model above.

# 1S0, 3S1, 3P0, 3P2
settings['GeV_to_MeV'] = np.array([1e-2,1e-2,1e-8,1e-8])

# Get constants from C++ module. Some constants are not acessable from the 
# python interface. To be able to reproduct results all constants are loaded
# and stored in the settings dictionary.
settings['c_Np']            = m.get_momentum_grid_points()
settings['c_mpi']           = m.get_mpi()
settings['c_Mn']            = m.get_Mn()
settings['c_Mp']            = m.get_Mp()
settings['c_gA']            = m.get_gA()
settings['c_fpi']           = m.get_fpi()

# Natural scale for lecs of certain dimensions

fac = (500.0/settings['lambda_cut'])
settings['GeV-2 scale'] = 0.13        # 10^4 GeV^{-2}
settings['GeV-4 scale'] = 0.63*fac**2 # 10^4 GeV^{-4}
settings['GeV-6 scale'] = 2.5*fac**4  # 10^4 GeV^{-6}   

# -------------------------------------
# Read experimental data
# -------------------------------------

# Read in the pandas dataframe
df = dm.dataframe_from_file(fname)  

# Get the chose observables
np = df.loc[df['TZ'] == 'NP'] # Take just np data
A = np.loc[np['Observable'] == 'DSG'] # Take just PB observable

# -------------------------------------
# Run the sampling code
# -------------------------------------
now = datetime.now()
dt_string = now.strftime("%H_%M_%S_%d-%m")
filename = '../output/sampl_' + dt_string + '.h5'

out = emcee_sampler(settings, data, filename, bayes.log_posterior)

# -------------------------------------
# Save output and settings
# -------------------------------------
d = {}
d['out']      = out
d['settings'] = settings
d['data']     = data
path_data     = '../output/sampl_'+ dt_string + '.npy'
np.save(path_data, d)
# To load the data: np.load(path_data,allow_pickle='TRUE').item()

# Print the runtime of the program
stop = tt.default_timer()
print("Program ended sucessfully. Runtime: " + str((stop-start)/60.0) + ' min')
