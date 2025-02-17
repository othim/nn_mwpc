import numpy as np
from scipy.special import genlaguerre
from scipy.special import gamma

# Load data
grid = np.loadtxt('grid.txt')
print(grid.shape)
p_grid = grid[:,0]
w_grid = grid[:,1]
pot = np.loadtxt('pot.txt').reshape(60,60)
print(pot.shape)


n  = 0
nn = 0
l  = 0
ll = 0

Omega = 24
mN    = 938.91/2 # This is wierd! Should be mN!?
b = np.sqrt(1/(mN*Omega))

def Rnl(p_grid,n,l):
    alpha = l+1.0/2.0
    
    x = p_grid**2*b**2
    L = genlaguerre(n,alpha)(x)

    return (-1)**n*np.sqrt(2*gamma(n+1)*b**3/(gamma(n+l+3.0/2.0)))*np.sqrt(x)**l*\
            np.exp(-x/2.0)*L

R_out = Rnl(p_grid,0,0)
R_in  = Rnl(p_grid,0,0)


me = 0
for i,p_out in  enumerate(p_grid):
    for j,p_in in enumerate(p_grid):
        me += w_grid[i]*w_grid[j]*p_in**2*p_out**2*R_out[i]*R_in[j]*pot[i,j]


print(f'me={me}')
