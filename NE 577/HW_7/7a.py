'''##### HW.7-a ##################################################'''
import numpy as np
import math
import matplotlib.pyplot as plot

#%% Constants 
nu = 1e-6   # kin' viscosity
mu = 1e-3   # viscosity
rho = 1e3   # density

#%% Grids
L_x = 0.04 # x
L_y = 0.02 # y

#%% Inputter
ny = 30         ## Input, no. of cells
nx = int(ny*2) 
h = L_x/nx      # resolution
 
#%% Droplet Data
rad = L_y/4
x_d = 0.02
y_d = 0.01

#%% Function Definition
def levelset (x,y):
    var = rad - math.sqrt( (x-x_d)**2 + (y - y_d)**2 )
    return var  

#%% Looper
x_n = np.linspace(0, L_x, nx+2)
y_n = np.linspace(0, L_y, ny)
LS_Field = np.zeros([len(x_n), len(y_n)])

x_mat = np.zeros([len(x_n),len(y_n)]) 
y_mat = np.zeros([len(x_n), len(y_n)])

for i in range(0,len(x_n)):
    for j in range(0,len(y_n)):
        y_mat[i,:] = y_n
        x_mat[:,j] = x_n
        LS_Field[i,j] = levelset(x_mat[i,j], y_mat[i,j]) 


#%% The Plotter
figure, axes = plot.subplots()
plot.contourf(x_mat, y_mat, LS_Field, 12, cmap ='Blues')
plot.title('Level Set Field') 
legend=plot.colorbar() 
legend.ax.set_title('Distance',fontsize = 10) 
plot.xlabel('Domain Length') 
plot.ylabel('Domain Height') 
interface = plot.Circle((x_d, y_d), rad, color ='blue', fill = False)
axes.add_artist(interface)
plot.show()

