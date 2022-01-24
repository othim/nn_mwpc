'''
    This program is an example of how to use the python interface to comute
    observables and phase shifts for some given potential models.

    - Imoprt the package and create and object with the potential you want

    - Print the LECs that the potential use, IT IS IN THIS ORDER YOU NEED
      TO ENTER THE VALUES IN THE LIST OF LECS GIVEN TO THE compute_observable
      FUNCTION!

    - To be sure that you have enetered them in the correct order, you can 
      check afterwards that the LEC values are as you expect.

    There are lots of more settings you can change, but that has to be done 
    in the C++ code and then recompile the package. 

    Oliver Thim 2021
'''


import nn_mwpc
import numpy as np 
import time
import matplotlib.pyplot as plt

# This function is not correct!! Not stable for all inputs!!

def blattToStapp(delta_minus_BB, delta_plus_BB, twoEpsilonJ_BB):                
    # Return phases in deg                                                      
    # Stapp convention (bar-phase shifts) in terms of Blatt-Biedenharn convention
    twoEpsilonJ = np.arcsin(np.sin(twoEpsilonJ_BB)*np.sin(delta_minus_BB - delta_plus_BB))  # mixing parameter
    delta_minus = 0.5*(delta_plus_BB + delta_minus_BB + np.arcsin(np.tan(twoEpsilonJ)/np.tan(twoEpsilonJ_BB)))
    delta_plus  = 0.5*(delta_plus_BB + delta_minus_BB - np.arcsin(np.tan(twoEpsilonJ)/np.tan(twoEpsilonJ_BB)))
    epsilon = 0.5*twoEpsilonJ                                                                                                                                    
    return delta_minus, delta_plus, epsilon 
# -----------
# Observables
# -----------

print("Constructing object")
#obj = nn_mwpc.nn_mwpc_interface("MWPC_LO_1",25,450.0,True,True)
obj = nn_mwpc.nn_mwpc_interface("WPC_LO",25,500.0,True,True)

print("Done")

print(obj.print_LECs_in_use())

# Just some test values
C1S0 = -0.112927/100.0
C3S1 = -0.087340/100.0

C1S0 = -0.1/100.0
gA2  = 1.29*1.29;
C3P0 = 0.1;
C3P2 = 0.1;


angles = np.linspace(1,179,50)
angles_l = angles.tolist()
energies = [10.0,50.0,200.0]

# Time the funtion call
start = time.time()
#obs_vector = obj.compute_observable("I 0000",angles_l,energies,[C1S0,C3P0,C3P2,C3S1,gA2])
obs_vector = obj.compute_observable_l("I 0000",angles_l,energies,[C1S0,C3S1,gA2])
end = time.time()

for E in np.linspace(1,30,30):

    phase_shift_vec = obj.compute_phase_shift(0,E,[C1S0,C3S1,gA2]);
    phase = phase_shift_vec[3]
    print(f'{E}   {phase*180/np.pi}')

print(obj.print_LEC_values()) # To confirm they are correct

#print(f'Total time: {1e3*(end-start):0.3f} ms,  Per energy: {1e3*(end-start)/len(energies):0.3f} ms')

print(obs_vector)

# ------------
# Phase shifts
# ------------

#model = "WPC_LO"
model = "MWPC_LO_1"
obj = nn_mwpc.nn_mwpc_interface(model,8,450.0,True,True)
print(obj.print_LECs_in_use())
T_lab = np.linspace(10,50,100) # MeV


C1S0 = 0
C3S1 = 0
gA2  = 0
C3P0 = 0
C3P2 = 0

if (model=="MWPC_LO_1"):
    # LECs from Jerry 
    C1S0 = -3.52087e-6
    C3S1 = -6.2791150e-6
    gA2  = 1.27*1.27;
    C3P0 = 4.78808625866e-11;
    C3P2 = -1.15148500523e-11;

    # Conversion to my conventions
    fac = (2/np.pi)*(2*np.pi)**3
    C1S0 = (C1S0 -3.70998077e-6)*(2/np.pi)*(2*np.pi)**3 # Seems correct
    C3S1 = C3S1*(2/np.pi)*(2*np.pi)**3 
    #gA2 = gA2*fac
    C3P0 = C3P0*fac
    C3P2 = C3P2*fac
elif(model=="WPC_LO"):
    C1S0 = -0.112927/100.0
    C3S1 = -0.087340/100.0
    gA2 = 1.289*1.289

# 1S0, 3S1, 3P0, 3P2
chn_i_list = [3]

i = 0
for chn_i in chn_i_list:

    phases_m = []
    phases_p = []
    phases_e = []
    for Tl in T_lab:
        phase_shift_vec = []
        if (model == "MWPC_LO_1"):
            phase_shift_vec = obj.compute_phase_shift(chn_i,Tl,[C1S0,C3P0,C3P2,C3S1,gA2]);
        elif(model =="WPC_LO"):
            phase_shift_vec = obj.compute_phase_shift(chn_i,Tl,[C1S0,C3S1,gA2]);
        #print(phase_shift_vec)

        if (chn_i == 3 or chn_i == 7):
            
            #a, b, c = blattToStapp(phase_shift_vec[1],phase_shift_vec[0],phase_shift_vec[2]*2.0);
            a, b, c = phase_shift_vec[1], phase_shift_vec[0], phase_shift_vec[2]
            phases_m.append(a)
            phases_p.append(b)
            phases_e.append(c)
        else:
            phases.append(phase_shift_vec[3])

    #print(f'Phase shifts: {phase_shift_vec}')
    #print(np.array(phases)*180/np.pi)
    
    fix, ax = plt.subplots(2,2)
    axs = ax.reshape(-1)
    axs[0].plot(T_lab,np.array(phases_m)*180/np.pi,'.-')
    axs[1].plot(T_lab,np.array(phases_p)*180/np.pi,'.-')
    axs[2].plot(T_lab,np.array(phases_e)*180/np.pi,'.-')
    axs[3].plot(T_lab, np.sin(np.array(phases_m)-np.array(phases_p)),'.-')

print(obj.print_LEC_values()) # To confirm they are correct
plt.show()
print('Exit OK')
