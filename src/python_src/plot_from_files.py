# This is a script that takes in files where these is 
# Two columsn of data and basically plots both
# of them in a graph to be able to compare



import numpy as np
import matplotlib.pyplot as plt
import subprocess

# Compile and run the ./obs file

#subprocess.call(["make obs"], shell=True, cwd="../cpp_src")
#subprocess.call('../cpp_src/obs')


# Read the data from the file names
data_1S0 = np.loadtxt('../../data/np_1S0_nijm1.txt')
data_3S1 = np.loadtxt('../../data/np_3S1_nijm.txt')

out_1S0 = np.loadtxt('../../data/out_1S0.txt')
out_3S1 = np.loadtxt('../../data/out_3S1.txt')

#print("Data")
#print(data_1S0)

#print("Output")
#print(out_1S0)

diff = data_1S0[:,1]-out_1S0[:,1];
energy = data_1S0[:,0]
print(diff)


# Plot the data
fig,ax = plt.subplots(1,1)
ax.plot(energy, diff)
#ax.plot(data_1S0[:,0], data_1S0[:,1], 'r')

# Plot the data
fig,ax = plt.subplots(1,1)
#ax.plot(energy, diff)
ax.plot(data_1S0[:,0], data_1S0[:,1], 'r')
ax.plot(out_1S0[:,0], out_1S0[:,1], 'g')


fig,ax = plt.subplots(1,1)
#ax.plot(energy, diff)
ax.plot(data_3S1[:,0], data_3S1[:,1], 'r')
ax.plot(out_3S1[:,0], out_3S1[:,2]%180, 'g')

plt.show()

## END ##
