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
import sys


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
    
def blattToStapp_corr(delta_m_BB,delta_p_BB,twoeps_BB):
    
    if (np.abs(delta_m_BB.imag) > 1e-12 or np.abs(delta_p_BB.imag) > 1e-12 or np.abs(twoeps_BB.imag) > 1e-12):
        print("Error,phases are complex")
        return np.NAN,np.NAN,np.NAN
    
    delta_m_BB = delta_m_BB.real
    delta_p_BB = delta_p_BB.real
    twoeps_BB  = twoeps_BB.real

    eps_BB = twoeps_BB/2.0

    cos2eps = np.cos(eps_BB)**2
    cos_2dp = np.cos(2*delta_p_BB)
    sin_2dp = np.sin(2*delta_p_BB)
    cos_2dm = np.cos(2*delta_m_BB)
    sin_2dm = np.sin(2*delta_m_BB)

    aR = cos2eps*cos_2dm + (1-cos2eps)*cos_2dp
    aI = cos2eps*sin_2dm + (1-cos2eps)*sin_2dp

    delta_m = 0.5*np.arctan2(aI,aR)
    
    aR = cos2eps*cos_2dp + (1-cos2eps)*cos_2dm
    aI = cos2eps*sin_2dp + (1-cos2eps)*sin_2dm
    delta_p = 0.5*np.arctan2(aI,aR)

    tmp  = 0.5*np.sin(2.0*eps_BB)
    aR   = tmp*(cos_2dm - cos_2dp)
    aI   = tmp*(sin_2dm - sin_2dp)
    tmp2 = (delta_p+delta_m)
    eps  = 0.5*np.arcsin(aI*np.cos(tmp2)-aR*np.sin(tmp2))

    return delta_m, delta_p, eps



def phase_shifts_Stapp(T11,T12,T22,p_on,mu,CORR_CONV):
    fac = 2*np.pi*1j*mu*p_on

    # Blatt-Biedenharn (BB) convention
    twoEpsilonJ_BB = np.arctan(2*T12/(T11-T22))	# mixing parameter
    delta_plus_BB  = -0.5*1j*np.log(1 - fac*(T11+T22)/2 + fac*(T12)/np.sin(twoEpsilonJ_BB))
    delta_minus_BB = -0.5*1j*np.log(1 - fac*(T11+T22)/2 - fac*(T12)/np.sin(twoEpsilonJ_BB))
 
    if (CORR_CONV):
        delta_minus, delta_plus, epsilon = blattToStapp_corr(delta_minus_BB, delta_plus_BB, twoEpsilonJ_BB)
    else:
        delta_minus, delta_plus, epsilon = blattToStapp(delta_minus_BB, delta_plus_BB, twoEpsilonJ_BB)

    return np.array([delta_minus, delta_plus, epsilon])

def cont_phase_shift_3S_D1(pm,pp,eps):
    # Make phase shift continous
    # Assume that the phase shift at infinity is zero
    # and that it varies continously

    n_pm = np.copy(pm)
    n_pp = np.copy(pp)
    n_eps = np.copy(eps)

    for i in range(1,len(pm)):
        # If we have a 180 degree discontinuity
        if (np.abs(pm[i]-pm[i-1])>150):
            # Add 180 degrees to all phases up until this point
            n_pm[0:i] = n_pm[0:i]+180

            # Change sign to the epsilon
            n_eps[0:i] = -n_eps[0:i]

    return n_pm, n_pp, n_eps

#
# This part of the code contains functions to compute phase shifts
# perturbatively in coupled channels in the Stapp convention.
#

#
# Define some functions
#

def f11(eps,d):
    return np.cos(2*eps)*np.exp(2*1j*d)

def f12(eps,d1,d2):
    return 1j*np.sin(2*eps)*np.exp(1j*(d1+d2))

#
# Derivatives f11
# 

# 1:st order

def deps_f11(eps,d):
    return -2*np.sin(2*eps)*np.exp(2*1j*d)

def dd_f11(eps,d):
    return 2*1j*np.cos(2*eps)*np.exp(2*1j*d)


# 2:nd order

def deps2_f11(eps,d):
    return -4*np.cos(2*eps)*np.exp(2*1j*d)

def dd_deps_f11(eps,d):
    return -4*1j*np.sin(2*eps)*np.exp(2*1j*d)

def dd2_f11(eps,d):
    return -4*np.cos(2*eps)*np.exp(2*1j*d)

# 3:rd order

def deps3_f11(eps,d):
    return 8*np.sin(2*eps)*np.exp(2*1j*d)

def dd2_deps_f11(eps,d):
    return 8*np.sin(2*eps)*np.exp(2*1j*d)

def dd_deps2_f11(eps,d):
    return -8*1j*np.cos(2*eps)*np.exp(2*1j*d)

def dd3_f11(eps,d):
    return -8*1j*np.cos(2*eps)*np.exp(2*1j*d)


#
# Derivatives f11
# 

# 1:st order

def deps_f12(eps,d1,d2):
    return 2*1j*np.cos(2*eps)*np.exp(1j*(d1+d2))

def dd_f12(eps,d1,d2):
    return -np.sin(2*eps)*np.exp(1j*(d1+d2))

# 2:nd order

def deps2_f12(eps,d1,d2):
    return -4*1j*np.sin(2*eps)*np.exp(1j*(d1+d2))

def dd_deps_f12(eps,d1,d2):
    return -2*np.cos(2*eps)*np.exp(1j*(d1+d2))

def dd2_f12(eps,d1,d2):
    return -1j*np.sin(2*eps)*np.exp(1j*(d1+d2))

# 3:rd order

def deps3_f12(eps,d1,d2):
    return -8*1j*np.cos(2*eps)*np.exp(1j*(d1+d2))

def dd2_deps_f12(eps,d1,d2):
    return -2*1j*np.cos(2*eps)*np.exp(1j*(d1+d2))

def dd_deps2_f12(eps,d1,d2):
    return 4*1j*np.sin(2*eps)*np.exp(1j*(d1+d2))

def dd3_f12(eps,d1,d2):
    return np.sin(2*eps)*np.exp(1j*(d1+d2))




#
# Define some more complicated functions that relates the 
# 


def g2_11(eps_0,d_0,eps_1,d_1):
    
    return (1/2)*(deps2_f11(eps_0,d_0)*eps_1**2 +  dd2_f11(eps_0,d_0)*d_1**2)+\
            dd_deps_f11(eps_1,d_0)*d_1*eps_1

def g3_11(eps_0,d_0,eps_1,d_1,eps_2,d_2):

    return deps2_f11(eps_0,d_0)*eps_1*eps_2 + dd_deps_f11(eps_0,d_0)*\
            (eps_1*d_2+eps_2*d_1) + dd2_f11(eps_0,d_0)*d_1*d_2 +\
            (1/6)*deps3_f11(eps_0,d_0)*eps_1**3 +\
            (1/2)*dd_deps2_f11(eps_0,d_0)*eps_1**2*d_1+\
            (1/2)*dd2_deps_f11(eps_0,d_0)*eps_1*d_1**2+\
            (1/6)*dd3_f11(eps_0,d_0)*d_1**3


def g2_12(eps_0,d1_0,d2_0,eps_1,d1_1,d2_1):

    return  (1/2)*deps2_f12(eps_0,d1_0,d2_0)*eps_1**2+\
            (1/2)*dd2_f12(eps_0,d1_0,d2_0)*d1_1**2+\
            (1/2)*dd2_f12(eps_0,d1_0,d2_0)*d2_1**2+\
            dd_deps_f12(eps_0,d1_0,d2_0)*eps_1*d1_1+\
            dd_deps_f12(eps_0,d1_0,d2_0)*eps_1*d2_1+\
            dd2_f12(eps_0,d1_0,d2_0)*d1_1*d2_1

def g3_12(eps_0,d1_0,d2_0,eps_1,d1_1,d2_1,eps_2,d1_2,d2_2):

    return  deps2_f12(eps_0,d1_0,d2_0)*eps_1*eps_2+\
            dd2_f12(eps_0,d1_0,d2_0)*d1_1*d1_2+\
            dd2_f12(eps_0,d1_0,d2_0)*d2_1*d2_2+\
            dd_deps_f12(eps_0,d1_0,d2_0)*(eps_1*d1_2+eps_2*d1_1)+\
            dd_deps_f12(eps_0,d1_0,d2_0)*(eps_1*d2_2+eps_2*d2_1)+\
            dd2_f12(eps_0,d1_0,d2_0)*(d1_1*d2_2+d1_2*d2_1)+\
            (1/2)*dd_deps2_f12(eps_0,d1_0,d2_0)*(eps_1**2*d1_1)+\
            (1/2)*dd_deps2_f12(eps_0,d1_0,d2_0)*(eps_1**2*d2_1)+\
            (1/2)*dd2_deps_f12(eps_0,d1_0,d2_0)*(eps_1*d1_1**2)+\
            (1/2)*dd3_f12(eps_0,d1_0,d2_0)*(d2_1*d1_1**2)+\
            (1/2)*dd2_deps_f12(eps_0,d1_0,d2_0)*(eps_1*d2_1**2)+\
            (1/2)*dd3_f12(eps_0,d1_0,d2_0)*(d1_1*d2_1**2)+\
            (1/6)*deps3_f12(eps_0,d1_0,d2_0)*eps_1**3+\
            (1/6)*dd3_f12(eps_0,d1_0,d2_0)*d1_1**3+\
            (1/6)*dd3_f12(eps_0,d1_0,d2_0)*d2_1**3
'''
def g3_12(eps_0,d1_0,d2_0,eps_1,d1_1,d2_1,eps_2,d1_2,d2_2):

    return  deps2_f12(eps_0,d1_0,d2_0)*eps_1*eps_2+\
            dd2_f12(eps_0,d1_0,d2_0)*d1_1*d1_2+\
            dd2_f12(eps_0,d1_0,d2_0)*d2_1*d2_2+\
            dd_deps_f12(eps_0,d1_0,d2_0)*(eps_1*d1_2+eps_2*d1_1)+\
            dd_deps_f12(eps_0,d1_0,d2_0)*(eps_1*d2_2+eps_2*d2_1)+\
            dd2_f12(eps_0,d1_0,d2_0)*(d1_1*d2_2+d1_2*d2_1)
'''   

def get_derivative_matrix(eps_0,d1_0,d2_0):
     return np.array([[deps_f11(eps_0,d1_0), dd_f11(eps_0,d1_0),0],\
            [deps_f12(eps_0,d1_0,d2_0),dd_f12(eps_0,d1_0,d2_0),dd_f12(eps_0,d1_0,d2_0)],\
            [deps_f11(eps_0,d2_0),0,dd_f11(eps_0,d2_0)]])


# These two functions are simplified, for the case where the NLO contribution
# vanishes

def phase_shift_pert_coup_N2_3LO_S(S11,S12,S22,eps_0,d1_0,d2_0):
    # Set up S
    S = np.array([S11,S12,S22])

    # Set up matrix
    A = get_derivative_matrix(eps_0,d1_0,d2_0)
    #print(S.shape)
    #print(A.shape)
    # Invert system to get phase shift corrections
    if np.linalg.cond(A) < 1/sys.float_info.epsilon:
        d = np.linalg.inv(A)@S
    else:
        print("Error, matrix singular. Returning 0")
        d = np.array([0,0,0])

    return np.array([d[1],d[2],d[0]])
    

def phase_shift_pert_coup_N2_3LO_T(T11,T12,T22,mu,p_on,eps_0,d1_0,d2_0):

    # Convert T to S
    fac = 2*np.pi*1j*mu*p_on
    S11 = -fac*T11
    S12 = -fac*T12
    S22 = -fac*T22
    # Call S
    return phase_shift_pert_coup_N2_3LO_S(S11,S12,S22,eps_0,d1_0,d2_0)


# These two are the full expressions

def phase_shift_pert_coup_N2LO_S(S11,S12,S22,eps_0,d1_0,d2_0,eps_1,d1_1,d2_1):
    # Set up S
    S = np.array([S11-g2_11(eps_0,d1_0,eps_1,d1_1),\
            S12-g2_12(eps_0,d1_0,d2_0,eps_1,d1_1,d2_1),\
            S22-g2_11(eps_0,d2_0,eps_1,d2_1)])

    # Set up matrix
    A = get_derivative_matrix(eps_0,d1_0,d2_0)
    #print(S.shape)
    #print(A.shape)
    # Invert system to get phase shift corrections
    if np.linalg.cond(A) < 1/sys.float_info.epsilon:
        d = np.linalg.inv(A)@S
    else:
        print("Error, matrix singular. Returning 0")
        d = np.array([0,0,0])

    return np.array([d[1],d[2],d[0]])


def phase_shift_pert_coup_N2LO_T(T11,T12,T22,mu,p_on,eps_0,d1_0,d2_0,eps_1,d1_1,d2_1):

    # Convert T to S
    fac = 2*np.pi*1j*mu*p_on
    S11 = -fac*T11
    S12 = -fac*T12
    S22 = -fac*T22
    # Call S
    return phase_shift_pert_coup_N2LO_S(S11,S12,S22,eps_0,d1_0,d2_0,eps_1,d1_1,d2_1)

def phase_shift_pert_coup_N3LO_S(S11,S12,S22,eps_0,d1_0,d2_0,eps_1,d1_1,d2_1,\
        eps_2,d1_2,d2_2):
    # Set up S
    g3 = g3_12(eps_0,d1_0,d2_0,eps_1,d1_1,d2_1,eps_2,d1_2,d2_2)
    
    S = np.array([S11-g3_11(eps_0,d1_0,eps_1,d1_1,eps_2,d1_2),\
            S12-g3,\
            S22-g3_11(eps_0,d2_0,eps_1,d2_1,eps_2,d2_2)])

    # Set up matrix
    A = get_derivative_matrix(eps_0,d1_0,d2_0)
    #print(S.shape)
    #print(A.shape)
    # Invert system to get phase shift corrections
    if np.linalg.cond(A) < 1/sys.float_info.epsilon:
        d = np.linalg.inv(A)@S
    else:
        print("Error, matrix singular. Returning 0")
        d = np.array([0,0,0])

    return np.array([d[1],d[2],d[0]])


def phase_shift_pert_coup_N3LO_T(T11,T12,T22,mu,p_on,eps_0,d1_0,d2_0,eps_1,d1_1,
        d2_1,eps_2,d1_2,d2_2):

    # Convert T to S
    fac = 2*np.pi*1j*mu*p_on
    S11 = -fac*T11
    S12 = -fac*T12
    S22 = -fac*T22
    # Call S
    return phase_shift_pert_coup_N3LO_S(S11,S12,S22,eps_0,d1_0,d2_0,eps_1,
            d1_1,d2_1,eps_2,d1_2,d2_2)
