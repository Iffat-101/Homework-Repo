

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm

Lx = 0.04
Ly = 0.02

rd = (1/2)*Ly/2

ny = 30
nx = 2*ny
dx = Lx/nx
dy = Ly/ny
h = dx
xd = 0.02
L_advect = 0.008

x = np.linspace(0+ dx/2, Lx- dx/2, nx) # Grid points in x-direction
y = np.linspace(-0.01+ dy/2, 0.01- dy/2, ny) # Grid points in y-direction
X, Y = np.meshgrid(x, y)

# Define time parameters
cfl_lim = 0.05
umax = 0.01875
dt= cfl_lim*(h/umax)
N_advect = int(L_advect/(umax*dt))
l2 = 1.0 # L2 norm initialization
n = 0 # iteration counter
ub = 0.0  # boundary velocity
epsilon = 1e-3
rhol = 1000
mul = 1e-3
rhog = 1
mug = 1e-5

phi = np.zeros((ny, nx+4))
sgnphi = np.zeros((ny,nx+4))
f  = np.zeros((ny, nx))
rho = np.zeros((ny, nx))
mu = np.zeros((ny, nx))
dpdx = np.zeros((ny, nx))
# nu = np.zeros((ny, nx))
M = 3

# Define velocity and pressure fields
u = np.zeros((ny+2, nx+2))  # u-velocity at cell edges in x-direction
v = np.zeros((ny+1, nx))  # v-velocity at cell edges in y-direction
p = np.zeros((ny, nx+2))  # pressure at cell centers
p_update = np.zeros((ny, nx+2)) 
ue =np.zeros((ny, nx))  # Exact velocity
A = u.copy() 
D = u.copy() 
u_star = u.copy()   # Intermediate u velocity in x-direction
v_star = v.copy()   # Intermediate v velocity in y-direction

# analytical solution
for j in range(0, ny):
    ue[j, :] = -187.5*y[j]**2 + 0.01875   

# initial condition
u[1:ny+1, 1:nx+1] = ue
u_star[1:ny+1, 1:nx+1] = ue 

            
# pressure initialization
p[:, :] = - 0.375*dx     

uc = np.zeros((ny, nx+2))
phi_star = np.zeros((ny, nx+4))
phi_update = np.zeros((ny, nx+4))
phi_phalf= np.zeros((ny, nx+4))
phi_mhalf= np.zeros((ny, nx+4))
Lphi= np.zeros((ny, nx+4))
Dxp = np.zeros((ny, nx))
Dxm = np.zeros((ny, nx))
Dxpp1 = np.zeros((ny, nx))
Dxmp1 = np.zeros((ny, nx))
Dxpm1 = np.zeros((ny, nx))
Dxmm1 = np.zeros((ny, nx))

 
# Level set distance function
for j in range(0, ny):
    for i in range(2, nx+2):
        phi[j, i] = rd -  np.sqrt((x[i-2]-xd)**2+(y[j])**2)

phi[:, 0] = phi[:, nx]  # Left wall
phi[:, 1] = phi[:, nx+1]  # Left wall
phi[:, nx+3] = phi[:, 2]  # Right wall
phi[:, nx+2] = phi[:, 3]  # Right wall

# Create the plot
plt.style.use("default")
plt.contourf(X, Y, phi[:, 2:nx+2], 500, cmap='magma')
plt.gca().set_aspect(1/1)
plt.xlabel('x')
plt.ylabel('y')
plt.title('2D contour of level set function')
plt.colorbar()
plt.savefig('line_contour_level_Set.png', bbox_inches = 'tight', dpi=100)
plt.show()

def property_update(phi_loc):
    # Calculate snoothed sign function
    for j in range(0, ny):
        for i in range(0, nx+4):
            if phi_loc[j,i] >= M*h:
               sgnphi[j,i] = 1
            elif phi_loc[j,i] <= -M*h:
               sgnphi[j,i] = -1
            else:
                sgnphi[j,i] = (phi_loc[j,i]/(M*h)) - ((1/np.pi)*np.sin((np.pi*phi_loc[j,i])/(M*h)))
                
    # Calculate the marlker function
    for j in range(0, ny):
        for i in range(0, nx):
            if phi_loc[j,i+2] < -M*h:
               f[j,i] = 0
            elif abs(phi_loc[j,i+2]) <= M*h:
               f[j,i]  = 0.5*(1+(phi_loc[j,i+2]/(M*h)) + ((1/np.pi)*np.sin((np.pi*phi_loc[j,i+2])/(M*h))))
            else:
                f[j,i]  = 1

    for j in range(0, ny):
        for i in range(0, nx):
            rho[j, i] = rhol*f[j, i] + rhog*(1 - f[j, i])
            
    for j in range(0, ny):
        for i in range(0, nx):
            mu[j, i] = mul*f[j, i] + mug*(1 - f[j, i])  
            
    for j in range(0, ny):
        for i in range(0, nx):
            rho[j, i] = rhol*f[j, i] + rhog*(1 - f[j, i])
            
    for j in range(0, ny):
        for i in range(0, nx):
            mu[j, i] = mul*f[j, i] + mug*(1 - f[j, i])
            
            
    nu = mu / rho   
    
    dpdx = -umax * mu * (2/(0.01)**2)
    
    return nu, dpdx, rho, mu, sgnphi

# Define functions for velocity and pressure updates
def inermediate_velocity_update(u, v, dt, dx, dy, rho, nu,p):
    un = u.copy()
    
    for j in range(1, ny+1):
        for i in range(1, nx+1):
            A[j, i] = (1/dx)*(((un[j, i+1] + un[j, i])/2)**2 -
                          ((un[j, i] + un[j, i-1])/2)**2)
            
    for j in range(1, ny+1):
        for i in range(1, nx+1):
            D[j, i] = (1/ (dx*dy))*(un[j+1, i] - 4*un[j, i] + un[j-1, i] + un[j, i+1]+ 
                          un[j, i-1])

    # Update intermediate velocities
    for j in range(1, ny+1):
        for i in range(1, nx+1):
           # u_star[j, i] = un[j, i] - dt *A[j, i]  + nu[j-1, i-1] * dt * D[j, i] - \
           #     dt * (1/rho[j-1, i-1]) * dpdx[j-1, i-1]
           
           u_star[j, i] = un[j, i] - dt *A[j, i]  + nu[j-1, i-1] * dt * D[j, i]
                         
                         
    return u_star

def pressure_update(u_star, v_star, dt, dx, dy, rho, epsilon,p):

    p_temp = p
    tolerance = 1
    
    while (tolerance > epsilon):
        for i in range(1,nx+1):
            p_update[0,i] = ((1/3)*(p_temp[0+1,i]+p_temp[0,i+1]+p_temp[0,i-1])) - ((dx/3)*(rho[0, i-1]/dt)*(u_star[1,i+1] - u_star[1,i]))
            p_update[ny-1,i] = ((1/3)*(p_temp[ny-1-1,i]+p_temp[ny-1,i+1]+p_temp[ny-1,i-1])) - \
                ((dx/3)*(rho[ny-1, i-1]/dt)*(u_star[ny,i+1] - u_star[ny,i])) 
        
        for j in range(1, ny-1):
            for i in range(1, nx+1):
                p_update[j, i] = ((1/4)*(p_temp[j+1,i]+p_temp[j-1,i]+p_temp[j,i+1]+p_temp[j,i-1])) - \
                ((dx/4)*(rho[j, i-1]/dt)*(u_star[j+1,i+1] - u_star[j+1,i]))
                        
        tolerance = (abs(p_temp[:,1:nx+1] - p_update[:,1:nx+1])).max()
        print('Tolerance: ',tolerance)
        
        p_update[:, nx+1] = p_update[:,2]
        p_update[:, 0] = p_update[:,nx-1]
 
        p_temp = p_update
        
    return p_update

def velocity_update(u_star, v_star, dt, dx, dy, rho, p_new):
    
    # Update velocities
    for j in range(1, ny+1):
        for i in range(1, nx+1):
            # Update u-velocity at cell edges in x-direction
            u[j, i] = u_star[j, i] - (1/rho[j-1, i-1])*(dt)*dpdx[j-1, i-1]
           
    return u

def switch(a, b):
    if abs(a) < abs(b):
        return a            
    else:
        return b
    
def Lphi_calculation(phi_loc, uc_loc):
    
    for j in range(0, ny):
        for i in range(0, nx):
            Dxp[j,i] = phi_loc[j,i+2] - phi_loc[j,i+1]
            Dxm[j,i] = phi_loc[j,i+1] - phi_loc[j,i]
            
    for j in range(0, ny):
        for i in range(0, nx):
            Dxpp1[j,i] = phi_loc[j,i+3] - phi_loc[j,i+2]
            Dxmp1[j,i] = phi_loc[j,i+2] - phi_loc[j,i+1]
            
    for j in range(0, ny):
        for i in range(2, nx+2):
            if 0.5*(uc_loc[j, i-1] + uc_loc[j, i-2]) > 0:
                  phi_phalf[j,i]= phi_loc[j,i] + 0.5*switch(Dxp[j,i-2], Dxm[j,i-2])
            else:
                  phi_phalf[j,i]= phi_loc[j,i+1] - 0.5*switch(Dxpp1[j,i-2], Dxmp1[j,i-2])
            
    for j in range(0, ny):
        for i in range(0, nx):
            Dxpm1[j,i] = phi_loc[j,i+2] - phi_loc[j,i+1]
            Dxmm1[j,i] = phi_loc[j,i+1] - phi_loc[j,i]
            
    for j in range(0, ny):
        for i in range(2, nx+2):
            if 0.5*(uc_loc[j, i-3] + uc_loc[j, i-2]) > 0:
                  phi_mhalf[j,i]= phi_loc[j,i-1] + 0.5*switch(Dxpm1[j,i-2], Dxmm1[j,i-2])
            else:
                  phi_mhalf[j,i]= phi_loc[j,i] - 0.5*switch(Dxp[j,i-2], Dxm[j,i-2])
    
        
    phi_phalf[:, 0] = phi_phalf[:, nx]  # Left wall
    phi_phalf[:, 1] = phi_phalf[:, nx+1]  # Left wall
    phi_phalf[:, nx+3] = phi_phalf[:, 2]  # Right wall
    phi_phalf[:, nx+2] = phi_phalf[:, 3]  # Right wall
    
    phi_mhalf[:, 0] = phi_mhalf[:, nx]  # Left wall
    phi_mhalf[:, 1] = phi_mhalf[:, nx+1]  # Left wall
    phi_mhalf[:, nx+3] = phi_mhalf[:, 2]  # Right wall
    phi_mhalf[:, nx+2] = phi_mhalf[:, 3]  # Right wall
    
    
    for j in range(0, ny):
        for i in range(2, nx+2):
            Lphi[j,i] = -uc_loc[j, i-2]*((phi_phalf[j,i] - phi_mhalf[j,i])/h)
    
    Lphi[:, 0] = Lphi[:, nx]  # Left wall
    Lphi[:, 1] = Lphi[:, nx+1]  # Left wall
    Lphi[:, nx+3] = Lphi[:, 2]  # Right wall
    Lphi[:, nx+2] = Lphi[:, 3]  # Right wall   
    
    return Lphi             

# Run simulation in iteration
while (n < N_advect):
    
    nu, dpdx, rho, mu, sgnphi = property_update(phi)
    
    # Solve Navier-Stokes equations
    u_star = inermediate_velocity_update(u, v, dt, dx, dy, rho, nu,p)
    
    u_star[:, 0] = u_star[:, nx-1]  # Left wall
    u_star[:, nx+1] = u_star[:, 2]  # Right wall
    u_star[0, :] = 2*ub - u_star[1, :]  # Bottom wall
    u_star[ny+1, :] = 2*ub - u_star[ny, :]  # Top wall
    
    p_new = pressure_update(u_star, v_star, dt, dx, dy, rho, epsilon,p)
    
    p_new[:, nx+1] = p_new[:,2]
    p_new[:, 0] = p_new[:,nx-1]
    # print(p)
    
    u_updated = velocity_update(u_star, v_star, dt, dx, dy, rho, p_new)
    
    # B.C. for updated velocity
    u_updated[:, 0] = u_updated[:, nx-1]  # Left wall
    u_updated[:, nx+1] = u_updated[:, 2]  # Right wall
    u_updated[0, :] = 2*ub - u_updated[1, :]  # Bottom wall
    u_updated[ny+1, :] = 2*ub - u_updated[ny, :]  # Top wall

    u = u_updated
    
    # L2 norm calculation
    l2 = (norm(abs(u[1:ny,0]-ue[1:ny,0])))/ny
    
    # B.C. for updated velocity
    u[:, 0] = u[:, nx-1]  # Left wall
    u[:, nx+1] = u[:, 2]  # Right wall
    u[0, :] = 2*ub - u[1, :]  # Bottom wall
    u[ny+1, :] = 2*ub - u[ny, :]  # Top wall
    
    for j in range(0, ny):
        for i in range(1, nx+1):
            uc[j,i] = 0.5*(u[j+1,i]+u[j+1,i+1])
        
    uc[:, 0] = uc[:, nx-1]  # Left wall
    uc[:, nx+1] = uc[:, 2]  # Right wall
    
    Lphi = Lphi_calculation(phi, uc)
    for j in range(0, ny):
        for i in range(2, nx+2):
            phi_star[j, i] = phi[j, i] + dt*Lphi[j,i]
            
    phi_star[:, 0] = phi_star[:, nx]  # Left wall
    phi_star[:, 1] = phi_star[:, nx+1]  # Left wall
    phi_star[:, nx+3] = phi_star[:, 2]  # Right wall
    phi_star[:, nx+2] = phi_star[:, 3]  # Right wall
    
    Lphi_update = Lphi_calculation(phi_star, uc)
    
    for j in range(0, ny):
        for i in range(2, nx+2):
            phi_update[j, i] = phi[j, i] + (dt/2)*(Lphi[j,i]+Lphi_update[j,i])
            
    phi_update[:, 0] = phi_update[:, nx]  # Left wall
    phi_update[:, 1] = phi_update[:, nx+1]  # Left wall
    phi_update[:, nx+3] = phi_update[:, 2]  # Right wall
    phi_update[:, nx+2] = phi_update[:, 3]  # Right wall
    
    phi = phi_update
    
    n = n+1

            

# Create the plot
plt.style.use("default")
plt.contourf(X, Y, phi_update[:, 2:nx+2], 500, cmap='rainbow')
plt.gca().set_aspect(1/1)
plt.xlabel('x')
plt.ylabel('y')
plt.title('2D contour of level set function_update')
plt.colorbar()
plt.savefig('line_contour_level_Set_update.png', bbox_inches = 'tight', dpi=100)
plt.show()

# Create filled contour plot using imshow
plt.imshow(u, extent=[0, nx*dx, 0, ny*dy], cmap='rainbow', origin='lower')
plt.colorbar()
plt.title('2D u-velocity Contour Plot')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('velocity_contour.png', bbox_inches = 'tight')
plt.show()

plt.plot(ue[0:ny, nx-1],y[0:ny],color='r', zorder=5, linewidth=1.5,label="True")
plt.title('Numerical Vs Analytical u-Velocity Plot')
plt.xlabel('x')
plt.ylabel('y')
plt.plot(u[1:ny+1, nx-1],y[0:ny],color='b', zorder=5, linewidth=1.5,label="True")
line_labels = ["Analytical", "Numerical"]
plt.figlegend( line_labels,  loc = 'center', borderaxespad=0.3, ncol=1, labelspacing=0., bbox_to_anchor=(0.78, 0.82), prop={'size': 10} )
plt.savefig(str(nx)+'n'+str(n)+'.png', bbox_inches = 'tight')
plt.show()
