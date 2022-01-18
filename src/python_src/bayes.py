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
        The logarithm of the posterior of the data given the parameters.
    """
    # The parameters (LECs) are here in the 10^4 GeV^{-2n} units. The units 
    # that are entered in the function to compute observables are MeV^{-2n}.

    return log_lik(param, settings, data) + log_prior(param, settings)

def log_lik_one_obs(q_on_shell, mpi, lambda_chi, c_bar, obs_e, obs_e_var, obs_th):
    """
        Returns the log likelihood of the given data point given the parameters.
        It is NOT normalized.
    """
    # The EFT expansion parameter
    Q = np.maximum(mpi,q_on_shell)/lambda_chi
    
    # The EFT model error from summing all contributions beyond
    # leading order.
    obs_th_var = (obs_e**2)*(c_bar**2)*(Q**4)/(1-Q**2)
    
    # The total variance from the model and theory error
    var = obs_e_var + obs_th_var

    # The UN-normalized log-contribution to the likelihood
    return -((obs_e - obs_th)**2)/(2*var2)
    


# Likelihood data given lecs, prior for lecs
def log_lik(param,settings,data):
    """
        This function returns the log likelihood of the data sigma_exp
        given the sigma_th computed with gien lecs.
    """
   
    # Convert the LECs to the correct units
    param_MeV = np.dot(param, settings['GeV_to_MeV'])


    energies = [] # List of energies
    mpi = settings['c_mpi']
    lambda_chi = settings['lambda_chi']
    c_bar = settings['c_bar']

    sum = 0 # All contributions from all data will be added

    # Loop over energies since the ALL observabels for a given energy can 
    # be computed efficiently since the saclay-amplitudes only needs to be
    # computed once.
    for Tlab in energies:
        
        # Compute saclay amplitudes for this energy
        sac_amp = ...
        q_on_shell = ...
        for obs in data['obs']: # Loop over all observables with same energy
            
            # Get the experimental and theoretical values of the observable in 
            # question
            obs_e = ...
            obs_e_var = ... 
            obs_th = obj.get_observable()
            
            # Compute the contribution to the likelihood
            sum -= log_lik_one_obs(q_on_shell, mpi, lambda_chi, c_bar, obs_e, \
                    obs_e_var, obs_th)
    
    return sum

def log_prior(param,settings):
    """
        Return the log prior pdf for the parameters.
    """
    # The prior for the lecs are a gausian prior
    # with width 5, that expect the lecs to be of natural size.
    sum = 0
    for i in range(len(param)):
        sum -= param[i]**2/ \
                (2*(settings['width_prior'] * settings['lecs_scale'][i])**2)
    return sum


def ppd():
    """
        This funstion computes the posteror predictive distribution.
    """


def bayes_lepage():
    """
        This function computes the lepage error-plot for a given observable
        from all the sampled parameters. A Bayesian Lepage plot!
    """

