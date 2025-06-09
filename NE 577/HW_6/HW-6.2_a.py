import numpy as np
import math
import matplotlib.pyplot as plot

'''Assignments'''
Lx = 0.015          # domain length
Nx = 80             # no. of elements
dns_l = 1000        # Liquid Density
dns_g = 1           # Vapor desnity

h = Lx/(Nx-1)          # element length
ln_ob = Lx/2           # liquid droplet 
r = Lx/4               # droplet radius
x = np.ones(Nx)
rho = np.ones(Nx)
phi = np.ones(Nx)

''' Defining the level-Set Function'''
def LS(x1):
    vr1 = min(abs(x1-(ln_ob-r)), abs(x1-(ln_ob +r)))
    if x1 > (ln_ob - r) and x1 < (ln_ob + r):      
            return vr1
    else:   
            return -vr1

''' Defining the Heaviside Function'''        
def HvSd(phi, M):
    vr2 = (1/2)*(1 + phi/(M*h) + (1/math.pi)*math.sin(math.pi *phi /(M*h)))
        # Defining the functions returner
    if phi< -1*M*h:
        return 0
    
    elif abs(phi) <= (1*M*h):
        return vr2
    
    elif phi > M*h:
        return 1

'''Numerical Mass'''
def drop1():              # for any element with phi > 0
    return dns_l*h

def drop2(rho):           # for any element
    return rho*h

def drop3(rho):           # for all the elements  
    weight = (rho - dns_g) /(dns_l - dns_g) # Defining weight
    return dns_l *weight *h

''' Looper & plot '''
for M in range (1,4):
    m1 = 0; m2 = 0; m3 = 0   # Initializing Interface Thickness 
    for i in range (0, len(x)):
        x[i] = 0.0+i*h
        phi[i] = LS(x[i])
        rho[i] = dns_l *HvSd(phi[i], M) + dns_g *(1 - HvSd(phi[i], M))
        m3 = drop3(rho[i]) + m3
        
        if phi[i] > 0:
           m1 = drop1() + m1
           m2 = drop2(rho[i]) + m2
    
    '''****Printer***'''
    print('{:4.3e}, {:4.3e}, {:4.3e}, {:1d}'.format(m1, m2, m3, M))
    plot.plot(x, rho, label='M={:1d}'.format(M))
    plot.yscale('log')                  # Defining the plotting type
    plot.xlabel('Domain length (m)')    # Defining the x-label
    plot.ylabel('\u03C1 (kg/m\u00b3)')  # Defining the y-label
                    # unicode character for rho
    plot.title('Density Distribution (kg/m\u00b3)'.format(r), fontsize = 10)
                                    # /u00b - unicode superscript
    plot.legend()   # Showing the legend
    # plot.grid()     # Showing the grid

plot.show()
       
   
        
    

