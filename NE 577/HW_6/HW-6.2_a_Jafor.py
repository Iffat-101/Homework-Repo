import numpy as np
import math
import matplotlib.pyplot as plt
L_x=0.015                         ## domain length
N_x=80                            ## no. of elements
rho_l= 1000
rho_g= 1
m_1= 1
m_2= 2
m_3= 3
h= L_x/(N_x-1)                     ## element length
obj_len=L_x/2                      ## object length is half of the domain length
r= L_x/4                           ## droplet radious
x= np.ones(N_x)
rho=np.ones(N_x)
phi=np.ones(N_x)

def level_set(x1):
    if x1 > (obj_len-r) and x1 < (obj_len+r):
            return min(abs(x1-(obj_len-r)), abs(x1-(obj_len +r)))
    else:
            return -1*min(abs(x1-(obj_len-r)), abs(x1-(obj_len+r)))
def heaviside(phi,M):
    if phi< -1*M*h:
        return 0
    elif abs(phi)<=(1*M*h):
        return (1/2)*(1+phi/(M*h)+(1/math.pi)*math.sin(math.pi*phi/(M*h)))
    elif phi>M*h:
        return 1    
def m_droplet_op1():              ## for option 1: for any element with phi>0
    return rho_l*h
def m_droplet_op2(rho):           ## option 2: for any element
    return rho*h
def m_droplet_op3(rho):           ## for all the elements  
    density_wt=(rho-rho_g)/(rho_l-rho_g) ## definition of density weight
    return rho_l*density_wt*h
for M in range (1,4):
    m1=0                          ## initialization of interface thickness m1,m2,m3
    m2=0
    m3=0
    for i in range (0, len(x)):
        x[i]=0.0+i*h
        phi[i]=level_set(x[i])
        rho[i]=rho_l*heaviside(phi[i], M)+rho_g*(1-heaviside(phi[i], M))
        m3=m_droplet_op3(rho[i])+m3
        if phi[i]>0:
           m1=m_droplet_op1()+m1
           m2=m_droplet_op2(rho[i])+m2
    print('{:4.3e}, {:4.3e}, {:4.3e}, {:1d}'.format(m1, m2, m3, M))
    plt.plot(x, rho,label='M={:1d}'.format(M))
    plt.yscale('log')
    plt.xlabel('domain length (m)')
    plt.ylabel('rho (kg/m3)')
    plt.legend()
    plt.title('value of density distribution (kg/m3), droplet radius= {:4.3e} m'.format(r), fontsize=8)
    plt.grid()

plt.show()
       
   
        
    

