import run_tests_on_module as tt
import sys
sys.path.append("/Users/toliver/Documents/phd") # Add path to library
import nn_mwpc
import numpy as np 
import matplotlib.pyplot as plt
plt.style.use('../plot_examples/plt_settings.mplstyle')




mesh_points = range(20)
mesh_points = np.array(mesh_points)*5 + 5
mesh_points = mesh_points.astype(int)
print(mesh_points)

max_inf_arr = []
max_fin_arr = []
tol = 1e-5
for pp in mesh_points:
    obj_inf = nn_mwpc.nn_mwpc_interface("MWPC_LO_J",2,500.0,6,False,True,True,pp,False)
    obj_fin = nn_mwpc.nn_mwpc_interface("MWPC_LO_J",2,500.0,6,False,True,True,pp,True)

    max_inf = tt.run_tests_phase(obj_inf,tol)
    max_fin = tt.run_tests_phase(obj_fin,tol)
    
    max_inf_arr.append(max_inf)
    max_fin_arr.append(max_fin)
    print(f' # points: {pp}, Max. rel. err. inf: {max_inf}, fin: {max_fin}')

for i,pp in enumerate(mesh_points):
    print(f' # points: {pp}, Max. rel. err. inf: {max_inf_arr[i]}, fin: {max_fin_arr[i]}')


fig,ax = plt.subplots(1,1)

ax.plot(mesh_points,max_inf_arr,linestyle='--',marker='o',label='inf mesh')
ax.plot(mesh_points,max_fin_arr,linestyle='--',marker='o',label='fin mesh')
ax.set_xlabel('\# mesh points')
ax.set_ylabel('abs. rel. err')
ax.set_yscale('log')
ax.legend()
plt.show()

