'''##### HW.7-c ##################################################'''
import numpy as np
import math
import matplotlib.pyplot as plt

#%% Defining Functions
def Analyt_u(y):
    var_analyt = -187.5 *y**2 + 0.01875
    return var_analyt

def LSet (x,y):
    var_Op_LS = rad - math.sqrt( (x -x_d)**2 + (y -y_d)**2 )  
    return var_Op_LS

def HSd (phi, M):
    if phi < -M*h:
        return 0
    elif abs(phi) <= M*h:
        var_HS = 0.5 * (1 + phi/(M*h) + (1/math.pi) *math.sin(math.pi *phi/(M*h))) 
        return  var_HS
    elif phi > M*h:
        return 1
    
def Op_L (u, phi_plus, phi_minus):
    var_mlt =  -u*(phi_plus - phi_minus)/(h)
    return var_mlt
 
#%% Constants
M = 3 
nu = 1e-6           # Kin' viscosity
mu_l = 1e-3         
mu_g = 1e-5         # Viscosities
rh_l = 1e3
rh_g = 1            # Densities
grd_P = -0.375      # Pressure Gradient, Analyt_u
v_mx = 0.01875      # Max' Velocity, Analyt_u

#%% Grids
L_1 = 0.04          # Domain Lenght
L_2 = 0.02          # Domain Width

#%% Inputter
ny = 50             ## Input, no. of cells
nx = int(ny*2) 
h = L_1/nx          # Resolution

#% The CFL
CFL = 0.999         ## Input CFL
dt = CFL*h /v_mx

#%%droplet info
rad = L_2/4
x_d = 0.02       # Droplet size
y_d = 0.01       # Droplet size @ y

#%% Variable Properties
phi     = np.zeros([nx+2, ny])
phi_str = np.zeros([nx+2, ny])
phi_n1  = np.zeros([nx+2, ny])

rho = np.zeros([nx+2, ny])
mu  = np.zeros([nx+2, ny])
x_n = np.linspace(0, L_1, nx+2)
y_n = np.linspace(0, L_2, ny)

x_mat = np.zeros([len(x_n), len(y_n)]) 
y_mat = np.zeros([len(x_n), len(y_n)])

D_plus = np.zeros([len(x_n), len(y_n)])
D_mins = np.zeros([len(x_n), len(y_n)])

#%% Analytical Velocity
u = np.zeros([len(x_n), len(y_n)])
Analyt_vel = np.linspace(-L_2/2, L_2/2, ny)
for i in range(0, len(x_n)):
    for j in range(0, len(y_n)):
        u[i,j] = Analyt_u(Analyt_vel[j])

#%% Advection Grid
for i in range(0, len(x_n)):
    for j in range(0, len(y_n)):
        y_mat[i,:] = y_n
        x_mat[:,j] = x_n
        phi[i,j] = LSet(x_mat[i,j], y_mat[i,j])
        
#%% Advector
for itr in range(0, 150):    # The printer loop
    
    # The Predictor    
    for i in range(0, len(x_n)):
        for j in range(0, len(y_n)):
            rho[i,j] = rh_l *HSd(phi[i,j], M) + rh_g *(1 - HSd(phi[i,j], M))
            mu[i,j]  = mu_l *HSd(phi[i,j], M) + mu_g *(1 - HSd(phi[i,j], M))
           
            if i < len(x_n)-1:
                phi_str[i,j] = phi[i,j] + dt * Op_L(u[i,j], phi[i+1,j],phi[i-1,j])
            else:
                phi_str[i,j] = phi[i,j] + dt * Op_L(u[i,j], phi[0,j],phi[i-1,j])

    # The Corrector
    for i in range(0, len(x_n)):
        for j in range(0,len(y_n)):
            phi_n1[i,j] = phi[i,j] + dt/2 * (Op_L(u[i,j],  phi[0,j], phi[i-1,j]) + Op_L(u[i,j],  phi_str[0,j], phi_str[i-1,j]))

    #%% The Plotter            
    figure, axes = plt.subplots()
    plt.contourf(x_mat, y_mat, phi, 10, cmap= 'Reds')
    plt.title('Level Set Field  \n Timestep = %i' % itr) 
    legend=plt.colorbar() 
    legend.ax.set_title('distance from interface',fontsize=8) 
    plt.xlabel('Domain Length') 
    plt.ylabel('Domain Height') 
    plt.show()
        
    phi = phi_n1        # Assigning new values
