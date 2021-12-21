#
#
# This file contains code that calculated the log likelihood, log prior
# and log posterior given some experimental data and some model for
# theoretically computing the observables.
#
#
# Oliver Thim, 2021 - 
import numpy as np
import time
from scipy.stats import norm

# Function to sample with emcee sampler
def log_posterior(param,settings,data):
    """
        This function returns the log prob. of the bayesian
        posterior distribution the lecs given the model settings
        in the dictionary settings.
        
        Input:
        param    - the parameters that are sampled over

        settings - settings for the run

        data     - contains the experimental data
        
        Output:
        The logarithm of the posterior of the data given the parameters
    """
    
    return log_likelihood(param, settings, data) +\ 
                log_prior(param, settings)

# Likelihood data given lecs, prior for lecs
def log_likelihood(Tlab_L,ko_L,sigma_th,sigma_exp,sigma_exp_err,c_bar,Lambda_chi):
    """
        This function returns the log likelihood of the data sigma_exp
        given the sigma_th computed with gien lecs.
    """

    
    # Old
    '''
    sum = 0
    for i,sigma in enumerate(sigma_th):
        Q = np.maximum(const.mpi,ko_L[i])/Lambda_chi

        model_err2 = (sigma_exp[i]**2)*(c_bar**2)*(Q**4)/(1-Q**2)
        var2 = sigma_exp_err[i]**2 + model_err2
        sum -= ((sigma - sigma_exp[i])**2)/(2*var2)
    '''
    return sum

def log_prior(param,settings):
    """
        Return the log prior pdf for the parameters.
    """
    # The prior for the lecs are a gausian prior
    # with width 5, that expect the lecs to be of natural size.
    sum = 0
    for i in range(len(param)):
        sum -= param[i]**2/(2*(settings['width_prior'] * settings['lecs_scale'][i])**2)
    return sum

