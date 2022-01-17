#
#
# This file contains code for performing sampling of a log posterior 
# distribution.
#
#
# Oliver Thim, 2021 -

import numpy as np
import emcee
import time
from multiprocessing import Pool
import multiprocessing


def emcee_sampler(settings,data,filename,log_posterior):
    """
        Function to that performs the emcee sampling for the desired
        parameters
        
        settings - the settings dictionary

        data     - the eperimental data included

        filename - the filename to save the data to

        Returns a list of the sampled lecs.
    """
    # Settings for sampler
    ndim = settings['s_ndim']
    nwalkers = settings['s_nwalkers']
    nsteps = settings['s_nsteps']
    initial = settings['s_initial'] # Initial position for walkers

    # Do the sampling in paralell using Pool
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)

    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, \
                log_posterior,args=[settings,data],pool=pool,backend=backend)
        
        # Run the steps
        pos, prob, state = sampler.run_mcmc(initial, nsteps, progress=settings['progress'])

    # Reshape samples
    samples = sampler.chain
    acc = sampler.acceptance_fraction
    # Return samples and current position
    return [samples, pos, acc]

def emcee_sampler2(settings,data,filename):
    """
        Function to that performs the emcee sampling for the desired
        parameters and also saves the auto correlation time during the
        sampling

        Returns a list of the sampled lecs.
    """
    # Settings for sampler
    ndim = settings['ndim']
    nwalkers = settings['nwalkers']
    nsteps = settings['nsteps']
    initial = data['initial'] # Initial guess for lecs

    # Do the sampling in paralell using Pool
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)

    autocorr = []

    from progress.bar import Bar

    class SlowBar(Bar):
        suffix = '%(remaining_hours)d hours remaining'
        @property
        def remaining_hours(self):
            return self.eta // 3600
        def set_prog(p):
            self.progress = p
    bar = SlowBar()

    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior,args=[settings,data],pool=pool,backend=backend)
        # Run the burn steps
        print("Run nsteps=" + str(nsteps))
        for sample in sampler.sample(initial, iterations=nsteps,progress=settings['progress']):
            # Only check convergence every 100 steps
            if sampler.iteration%10:
                continue
            bar.set_prog(7)
            bar.next()
            # Compute the autocorrelation time so far
            # Using tol=0 means that we'll always get an estimate even
            # if it isn't trustworthy
            tau = sampler.get_autocorr_time(tol=0)
            autocorr.append(np.mean(tau))
            print(str(np.max(tau)) + ' ' + str(sampler.iteration))
            # Check convergence
            converged = np.all(tau * 100 < sampler.iteration)


            if converged:
                print('converged')
                break
    # Reshape samples
    bar.finish()
    samples = sampler.chain.reshape((-1, ndim))
    pos = 0
    # Return samples and current position
    return [samples, pos, autocorr]

# Function to produce predictions from the distributions of the inferred parameters

def loop_helper(sample):
    lecs = {}
    lecs['C1S0'] = sample[0]*settings_1['lecs_scale'][0]
    lecs['C3S1'] = sample[1]*settings_1['lecs_scale'][1]
    lecs['C3P0'] = sample[2]*settings_1['lecs_scale'][2]
    lecs['C3P2'] = sample[3]*settings_1['lecs_scale'][3]

    # compute all cross sections.
    sigma_l = []
    for i,Tlab in enumerate(T_lab_L_1):
        S_1, ko = cs.compute_S_1_new(Tlab,data_1['NN_channels'],settings_1['Lambda'],lecs,settings_1['Np'] ,data_1['pm'],data_1['w'])
        sigma = cs.compute_cross_section(data_1['NN_channels'], S_1, ko)

        if PRED_ERR_1:
            # Add model error
            Q = np.maximum(const.mpi,ko)/settings_1['Lambda_chi']
            model_err2 = (data_1['sigma_exp'][i]**2)*(c_bar**2)*(Q**4)/(1-Q**2)
            r = norm.rvs(size=1,scale=np.sqrt(model_err2))[0]
            sigma = sigma + r

        sigma_l.append(sigma)
    return sigma_l

def compute_predictions(settings, data, rand_samples, T_lab_L,PAR=True,PRED_ERR = True, number_of_samples=1):
    """
        Compute predictions for phase shifts using samples lecs
        in samples.

        Return an array of length (len(samples), len(T_lab_L)) of
        computed cross sections.
    """
    print('Number of samples: ',len(rand_samples)*number_of_samples)
    global settings_1
    settings_1 = settings
    global data_1
    data_1 = data
    global T_lab_L_1
    T_lab_L_1 = T_lab_L
    global PRED_ERR_1
    PRED_ERR_1 = PRED_ERR

    global c_bar
    c_bar      = settings['c_bar']
    sigma_pred = []

    # Make number_of_samples copys of rand_samples and
    # stack them on each other.
    if PRED_ERR:
        rand_samples = np.tile(rand_samples,(number_of_samples,1))

    if PAR:
        a_pool = multiprocessing.Pool()

        result = a_pool.map(loop_helper, rand_samples)
        return np.array(result)
    else:
        for s,sample in enumerate(rand_samples):

            if not(s%100):
                print('On sample nr: ', s)
            lecs = {}
            lecs['C1S0'] = sample[0]*settings['lecs_scale'][0]
            lecs['C3S1'] = sample[1]*settings['lecs_scale'][1]
            lecs['C3P0'] = sample[2]*settings['lecs_scale'][2]
            lecs['C3P2'] = sample[3]*settings['lecs_scale'][3]

            # compute all cross sections.
            sigma_l = []
            for i,Tlab in enumerate(T_lab_L):
                S_1, ko = cs.compute_S_1_new(Tlab,data['NN_channels'],settings['Lambda'],lecs,settings['Np'] ,data['pm'],data['w'])
                sigma = cs.compute_cross_section(data['NN_channels'], S_1, ko)

                if PRED_ERR:
                    # Add model error
                    Q = np.maximum(const.mpi,ko)/settings['Lambda_chi']
                    model_err2 = (data['sigma_exp'][i]**2)*(c_bar**2)*(Q**4)/(1-Q**2)
                    r = norm.rvs(size=1,scale=np.sqrt(model_err2))[0]
                    sigma = sigma + r
                sigma_l.append(sigma)

            sigma_pred.append(sigma_l)
        return np.array(sigma_pred)
