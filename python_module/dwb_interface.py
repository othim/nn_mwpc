'''
    This is a file that contains routines for the python interface for the
    nn_mwpc code to easily

    - Test functionality

    - Compute and plotobservables 

    - Compute and plot phase shifts

    - Compute and plot bound state energies

    Oliver Thim, 2023
'''
import numpy as np
import matplotlib.pyplot as plt
import time
import nn_mwpc



def get_pwa93_phase_shifts(chn):
    path = '/Users/toliver/Documents/phd/projects/lo_lepage/data/pwa93/'
    filename = 'np_' + chn + '_pwa93.txt' # chn is the LS term
    
    try:
        data = np.loadtxt(path + filename)
    except:
        print('No pwa phase shifts for this channel')
        data = []
    
    return data

def plot_pwa93_phase_shifts(chn,ax):

    try:
        data = get_pwa93_phase_shifts(chn)

        x_axis = data[0:-1:2,0] 
    
        idx = [1,2,3]

        if (data.shape[1]>2): # coupled channel 
            ax = ax.reshape(-1)        
            for i in range(3):
            
                if (i==0):
                    ax[i].plot(x_axis, data[0:-1:2,idx[i]],color='black',\
                            linestyle='--',label='pwa93')
            
                ax[i].plot(x_axis, data[0:-1:2,idx[i]],color='black',\
                        linestyle='--')
        else:
            ax.plot(x_axis, data[0:-1:2,1],color='black',\
                    linestyle='--',label='pwa93')
    except:
        print('No pwa phase shifts for this channel')


def plot_phase_shifts(Tlab,phase_shifts, chn_index, coupled, LS_term, DATA_DIR):
    '''
        

    '''

    creal = 'C0'
    cimag = 'C1'
    # print phase shifts to file or to the terminal
    if (coupled):
        
        fig, ax = plt.subplots(2,2,figsize=(9,7))
        
        plot_pwa93_phase_shifts(LS_term,ax)

        ax[0,0].plot(Tlab, np.real(phase_shifts[:,0]),color=creal,label='real')
        ax[0,0].plot(Tlab, np.imag(phase_shifts[:,0]),color=cimag,label='imag')
        
        ax[0,0].set_ylabel(r'$\delta_m$ (deg)')

        ax[0,1].plot(Tlab, np.real(phase_shifts[:,1]),color=creal)
        ax[0,1].plot(Tlab, np.imag(phase_shifts[:,1]),color=cimag)
        
        ax[0,1].set_ylabel(r'$\delta_p$ (deg)')
        
        ax[1,0].plot(Tlab, np.real(phase_shifts[:,2]),color=creal)
        ax[1,0].plot(Tlab, np.imag(phase_shifts[:,2]),color=cimag)
        
        ax[1,0].set_ylabel(r'$\epsilon$ (deg)')
        
        ax[1,0].set_xlabel(r'$T_\mathrm{lab}$ (MeV)')
        ax[0,1].set_xlabel(r'$T_\mathrm{lab}$ (MeV)')
        
        ax[1,1].set_visible(False)
        
        fig.suptitle(LS_term)
        fig.tight_layout()
        return fig,ax
    else:
        fig, ax = plt.subplots(1,1,figsize=(6.4, 4.8))

        plot_pwa93_phase_shifts(LS_term,ax)
        
        ax.plot(Tlab, np.real(phase_shifts[:,3]),color=creal,label='real')
        ax.plot(Tlab, np.imag(phase_shifts[:,3]),color=cimag,label='imag')
        
        ax.set_xlabel(r'$T_\mathrm{lab}$ (MeV)')
        ax.set_ylabel(r'$\delta$ (deg)')

        fig.suptitle(LS_term)
        fig.tight_layout()
        
        return fig,ax




def phase_shift(T,p_on,mu):
    
    return -0.5*1j*np.log(1-2*np.pi*1j*mu*p_on*T)

def d_phase_shift(T,p_on,mu):

    return -np.pi*mu*p_on/(1.0-2.0*1j*np.pi*mu*p_on*T)

def d2_phase_shift(T,p_on,mu):

    return -2.0*1j*(np.pi*mu*p_on)**2/(1.0-2.0*1j*np.pi*mu*p_on*T)**2

def d3_phase_shift(T,p_on,mu):
    
    return -(2*np.pi*mu*p_on)**3/(1.0-2.0*1j*np.pi*mu*p_on*T)**3

def phase_shift_pert_NLO(T0,T1,p_on,mu):
    
    
    d1 = d_phase_shift(T0,p_on,mu)*T1
    #print(d_phase_shift(T0,p_on,mu))
    #print(T1)
    #print(d1)
    #print("\n\n")
    return d1

def phase_shift_pert_N2LO(T0,T1,T2,p_on,mu):
    
    tmp1 = d_phase_shift(T0,p_on,mu)*T2 
    tmp2 = d2_phase_shift(T0,p_on,mu)*T1**2/2.0
    d2 = tmp1 + tmp2

    return d2

#def phase_shift_pert_N2LO(T0,T1,T2,p_on,mu):
#    
#    tmp1 = d_phase_shift(T0,p_on,mu)*T2 
#    d1 = phase_shift_pert_NLO(T0,T1,p_on,mu)
#    d2 = tmp1 - 1j*d1**2
#
#    return d2
    
def phase_shift_pert_N3LO(T0,T1,T2,T3,p_on,mu):
    
    tmp1 = d_phase_shift(T0,p_on,mu)*T3
    tmp2 = d2_phase_shift(T0,p_on,mu)*T1*T2
    tmp3 = d3_phase_shift(T0,p_on,mu)*T1**3/6

    return tmp1 + tmp2 + tmp3


def blattToStapp(delta_minus_BB, delta_plus_BB, twoEpsilonJ_BB):

    # Stapp convention (bar-phase shifts) in terms of Blatt-Biedenharn convention
    twoEpsilonJ = np.arcsin(np.sin(twoEpsilonJ_BB)*np.sin(delta_minus_BB - delta_plus_BB))	# mixing parameter
    delta_minus = 0.5*(delta_plus_BB + delta_minus_BB + np.arcsin(np.tan(twoEpsilonJ)/np.tan(twoEpsilonJ_BB)))
    delta_plus  = 0.5*(delta_plus_BB + delta_minus_BB - np.arcsin(np.tan(twoEpsilonJ)/np.tan(twoEpsilonJ_BB)))
    epsilon     = 0.5*twoEpsilonJ

    return delta_minus, delta_plus, epsilon
    
def phase_shifts_Stapp(T11,T12,T22,p_on,mu):
    fac = 2*np.pi*1j*mu*p_on

    # Blatt-Biedenharn (BB) convention
    twoEpsilonJ_BB = np.arctan(2*T12/(T11-T22))	# mixing parameter
    delta_plus_BB  = -0.5*1j*np.log(1 - fac*(T11+T22)/2 + fac*(T12)/np.sin(twoEpsilonJ_BB))
    delta_minus_BB = -0.5*1j*np.log(1 - fac*(T11+T22)/2 - fac*(T12)/np.sin(twoEpsilonJ_BB))
 
    
    delta_minus, delta_plus, epsilon = blattToStapp(delta_minus_BB, delta_plus_BB, twoEpsilonJ_BB)
    
    return np.array([delta_minus, delta_plus, epsilon])



