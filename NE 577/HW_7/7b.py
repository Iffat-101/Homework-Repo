'''##### HW.7-b ##################################################'''
import numpy as np
import math
import matplotlib.pyplot as plot
   
#%% Function Definitions
def LSt (x,y):
    var_LS = radius - math.sqrt( (x-x_d)**2 + (y - y_d)**2 ) 
    return var_LS

def HSd (phi, M):
    if phi< -M*h:
        return 0
    elif abs(phi) <= M*h:
        var_HS = 0.5 * (1+phi/(M*h) + (1/math.pi)*math.sin(math.pi*phi/(M*h)))
        return   var_HS
    elif phi > M*h:
        return 1

#%% Defining Constants 
nu = 1e-6       # Kin' viscosity
mu_l = 1e-3     # Viscosity liquid
mu_g = 1e-5     # Viscosity gas
rh_l = 1e3      # Density
rh_g = 0.6
M = 3
grd_P = -0.375

#%% Grids
L_x = 0.04 # x
L_y = 0.02 # y

#%% Inputter
ny = 50         ## Input, no. of cells
nx = int(ny*2) 
h = L_x/nx      # Resolution

#%% About Droplet
radius = L_y/4
x_d = 0.02
y_d = 0.01

##### Property Array
phi = np.zeros([nx+2,ny])
rho = np.zeros([nx+2,ny])
mu = np.zeros([nx+2,ny])

x_n = np.linspace(0,L_x,nx+2)
y_n = np.linspace(0,L_y, ny)  # Taking for boundary  nodes
level_set_field = np.zeros([len(x_n), len(y_n)])

x_mat = np.zeros([len(x_n),len(y_n)]) 
y_mat = np.zeros([len(x_n), len(y_n)])

#%% The Looper
for i in range(0,len(x_n)):
    for j in range(0,len(y_n)):
        y_mat[i,:] = y_n
        x_mat[:,j] = x_n
        phi[i,j] = LSt(x_mat[i,j], y_mat[i,j]) 
        rho[i,j] = rh_l*HSd(phi[i,j], M) + rh_g*(1-HSd(phi[i,j], M))
        mu[i,j]  = mu_l*HSd(phi[i,j], M) + mu_g*(1-HSd(phi[i,j], M))

#%% The Plotter
## Density
plot.contourf(x_mat, y_mat, rho, cmap = 'Blues')
legend = plot.colorbar()
plot.title('Density Distribution')
legend.ax.set_title(r' $\frac{kg}{m^3}  $') # the r' means latex format
plot.xlabel('Domain Length') 
plot.ylabel('Domain Height')
plot.show()

## Viscosity
plot.contourf(x_mat, y_mat, mu, cmap = 'Greens')
plot.title('Viscosity Disstribution')
legend = plot.colorbar()
legend.ax.set_title('$ Pa \cdot s $')
plot.xlabel('Domain Length') 
plot.ylabel('Domain Height')
plot.legend
plot.show()
