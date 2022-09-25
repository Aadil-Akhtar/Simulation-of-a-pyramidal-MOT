from re import U, X
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import random
import pandas as pd
import seaborn as sns
import math
import scipy
from scipy.stats import norm

# Importing all the necessary libraries

# number of atoms
num = 100

# the constants of the system

g_ = 3/2 
mu_B = 9.2740100783e-24 
h_bar = 1.055e-34 # reduced Planck's constant
gamma = 0.08 
det = (-2*22/7)*(110e3)
pi = 22.0/7.0     # pi
M = 133*1.67e-27  # mass of the atom
T = 300           # temperature in K
k = 1.38e-23      # Boltzmann constant
g = 9.8           # gravity
Kwave = (2*22)/(7*852e-9) # what is this?
del_t= 0.1/(2*pi*5.332e6) # time step
mom = (h_bar*Kwave/M) # momentum of the wave

Gamma = 61.542e6 #?
W1 = 1.0 #?
W2 = 1.0 #?
S = 9.0 #?
# A = 0.00008
A = 0.8 #?
det = (-2*22/7)*(4e3) # why 2 times det?
omega = 1*A*g_*mu_B*gamma/h_bar # the zeeman splitting factor



def spherical():
  u = np.random.random_sample()*2. - 1.0  # uniform random number between -1 and 1
  x = np.random.random_sample() - 0.5
  y = np.random.random_sample() - 0.5
  z = np.random.random_sample() - 0.5
  mag = np.sqrt(x**2 + y**2 + z**2)
  
  x=x/mag
  y=y/mag
  z=z/mag
  
  c = cbrt(u)  # for making almost constant volume denisty
  
  c = c*(1e-3)  # scaling to our choice 
  
  return x*c, y*c, z*c


# initializing the positions of the atoms
x_lst, y_lst, z_lst = [], [], []
for i in range(num):
  x, y, z = spherical()
  x_lst.append(5*x)
  y_lst.append(5*y)
  z_lst.append(5*z)


# plotting the points in space

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plotting the points

x_ell = 5e-3*np.outer(np.cos(np.linspace(0, 2 * np.pi, 100)), np.sin(np.linspace(0, np.pi, 100)))  # ellipsoid because of the density
y_ell = 5e-3*np.outer(np.sin(np.linspace(0, 2 * np.pi, 100)), np.sin(np.linspace(0, np.pi, 100)))  # ellipsoid because of the density
z_ell = (5e-3)*np.outer(np.ones(np.size(np.linspace(0, 2 * np.pi, 100))), np.cos(np.linspace(0, np.pi, 100))) # ellipsoid in z direction different because of the symmetry ??? 
"""
ax = plt.axes(projection='3d')  # 3d plot
ax.scatter3D(x_lst, y_lst, z_lst, s=2)  # s is the size of the points
ax.plot_surface(x_ell, y_ell, z_ell, alpha=0.2)  # alpha is the transparency

"""

#plt.show()
"""

"""




# Uptil here, we have plotted the points in space.

# clearing the figureuptil now  using plt.clf() 


# Now, we will generate random numbers for the velocities of the atoms.



#standard deviation of the velocities
sigma = np.sqrt( (k*T)/M ) # this is important because we want to have a gaussian distribution???



from scipy.stats import maxwell # for generating random numbers from a Maxwell distribution

u = maxwell.rvs(loc=0, scale=sigma, size=num) 

   
count, bins, ignored = plt.hist(u, 100, density=True, label='histogram')  # plotting the histogram of the random numbers  
plt.legend(loc='upper right') # legend
# plt.title('Initial speed distribution')
# plt.savefig("Initial speed distribution", dpi=1000)
# plt.show()



theta = np.pi*np.random.rand(num) 
phi = 2*np.pi*np.random.rand(num) 

ux = u*np.sin(theta)*np.cos(phi)  # x component of the velocity
uy = u*np.sin(theta)*np.sin(phi)  # y component of the velocity
uz = u*np.cos(theta)              # z component of the velocity

# ux,uy,uz are arrays/lists of the velocities of the atoms.


#Now, we will let the system evolve with time.

time_itt = 10000 # number of time steps
#Total time
time_tot = del_t*time_itt 
# print(time_tot) 

# initializing the lists for the positions and the velocities
x_pos, x_sp = [x_lst], [ux] 
y_pos, y_sp = [y_lst], [uy]
z_pos, z_sp = [z_lst], [uz]

#print(ux)
#print(x_sp[0][3])


for i in range(time_itt):
    dum_spx, dum_posx, dum_spy, dum_posy, dum_spz, dum_posz = [], [], [], [], [], []

    
    for j in range(num):
        
        x_sp[i][j] #
        p1 = (Gamma*W1*S/2)/(1 + W1*S + (4/Gamma**2)*(det - Kwave*x_sp[i][j] - omega*x_pos[i][j])**2 )
        p2 = (Gamma*W1*S/2)/(1 + W1*S + (4/Gamma**2)*(det + Kwave*x_sp[i][j] + omega*x_pos[i][j])**2 )

        dum_spx.append(x_sp[i][j] + mom*(p1 - p2)*del_t)
        dum_posx.append(x_pos[i][j] + x_sp[i][j]*del_t + mom*(p1 - p2)*(0.5*del_t**2))
        
        
       
          
          
            
        p3 = (Gamma*W1*S/2)/(1 + W1*S + (4/Gamma**2)*(det - Kwave*y_sp[i][j] - omega*y_pos[i][j])**2 )
        p4 = (Gamma*W1*S/2)/(1 + W1*S + (4/Gamma**2)*(det + Kwave*y_sp[i][j] + omega*y_pos[i][j])**2 )    

        dum_spy.append(y_sp[i][j] + mom*(p3 - p4)*del_t)
        dum_posy.append(y_pos[i][j] + y_sp[i][j]*del_t + mom*(p3 - p4)*(0.5*del_t**2))


        
        
                        

        p5 = (Gamma*W1*S/2)/(1 + W1*S + (4/Gamma**2)*(det - Kwave*z_sp[i][j] - 2*omega*z_pos[i][j])**2 )
        p6 = (Gamma*W1*S/2)/(1 + W1*S + (4/Gamma**2)*(det + Kwave*z_sp[i][j] + 2*omega*z_pos[i][j])**2 )
        dum_spz.append(z_sp[i][j] + mom*(p5 - p6)*del_t + g*del_t)
        dum_posz.append(z_pos[i][j] + z_sp[i][j]*del_t + mom*(p5 - p6)*(0.5*del_t**2) + 0.5*g*del_t**2)            
            

    
    
    
    x_sp.append(dum_spx)
    x_pos.append(dum_posx)
    y_sp.append(dum_spy)
    y_pos.append(dum_posy)
    z_sp.append(dum_spz)
    z_pos.append(dum_posz)


# why are all atoms absorbing the photon???
  

plt.scatter(x_lst, z_lst, s=4, c="b",)
plt.scatter(x_pos[-1], z_pos[-1], s=2, c='r',)
ax = plt.gca()
# plt.xlim([-5e-3, 5e-3])
# plt.ylim([-5e-3, 5e-3])

plt.plot(5e-3*np.cos(np.linspace(0, 2*np.pi, 100)) , (5e-3)*np.sin(np.linspace(0, 2*np.pi, 100)), linestyle='--')
ax.set_aspect('auto')
# plt.title('Initial position(blue) vs Final position(red) XZ')
plt.savefig("Initial position(blue) vs Final position(red) XZ", dpi=1000)
#plt.show()
plt.clf()


v = [0]*num
for q in range(num):
    v[q] = np.sqrt((x_sp[time_itt-1][q])**2 + (y_sp[time_itt-1][q])**2 + (z_sp[time_itt-1][q])**2)
    

# r = [0]*num

# for w in range(num):
#     r[w] = np.sqrt((np.array(x_pos[-1][w]))**2 + (np.array(y_pos[-1][w]))**2 + (np.array(z_pos[-1][w]))**2)   


# print(min(v))
# print(min(r))

# print(v.index(min(v)))
# print(r.index(min(r)))



count, bins, ignored = plt.hist(v, 100, density=True, label='histogram')
plt.legend(loc='upper right')
# plt.xlim
plt.title('Final speed distribution')
plt.savefig("Final speed distribution", dpi=1000)
# plt.show()
# plt.clf()
#plt.show()

count, bins, ignored = plt.hist(u, 100, density=True, label='histogram')
plt.legend(loc='upper right')
plt.title('Initial speed distribution')
plt.savefig("Initial speed distribution", dpi=1000)

# plt.show()

# plt.clf()

# avg = np.mean(v)
# var = np.var(v)
# sd = np.sqrt(var)

plt.subplot(1,2,1)

h =plt.hist2d(x_lst, z_lst, 20)
ax = plt.gca()
# ax.set_aspect(1)
plt.colorbar(h[3])
# plt.xticks(fontsize=15, rotation=90)
# plt.yticks(fontsize=15, rotation=0)
# plt.xlim([-0.0105, 0.0105])
# plt.ylim([-0.01, 0.01])
plt.xlim([-0.0105, 0.0105])
plt.ylim([-0.01, 0.01])
plt.title('Initial Position 2D Histogram')
plt.savefig("Initial Position 2D Histogram", dpi=1000)


#plt.clf()
plt.subplot(1,2,2)
h =plt.hist2d(x_pos[-1], z_pos[-1], 20)
ax = plt.gca()
# ax.set_aspect(1)
plt.colorbar(h[3])
plt.xticks(rotation=90)
# plt.xlim([-5e-3, 5e-3])
# plt.ylim([-5e-3, 5e-3])
plt.xlim([-0.0105, 0.0105])
plt.ylim([-0.01, 0.01])
plt.title('Final Position 2D Histogram')
plt.savefig("Final Position 2D Histogram", dpi=1000)
plt.show()



#plt.scatter(x_pos[-1], z_pos[-1], s=4, c='b') # plotting the points in space
#ax = plt.gca() # getting the current axis
#plt.plot(5e-3*np.cos(np.linspace(0, 2*np.pi, 100)) , (5e-3)*np.sin(np.linspace(0, 2*np.pi, 100)), linestyle='--') # cross -section ??? 
#ax.set_aspect('auto') # setting the aspect ratio
#plt.show()

#ax = plt.axes(projection='3d')  # 3d plot
#ax = plt.axes(projection='3d')  # 3d plot
#ax.scatter3D(x_lst, y_lst, z_lst, s=2, c='b')  # s is the size of the points

#ax.plot_surface(x_ell, y_ell, z_ell, alpha=0.2)  # alpha is the transparency


#ax.scatter3D(x_pos[-1], y_pos[-1], z_pos[-1], s=2, c='r')  # s is the size of the points
#ax.plot_surface(x_ell, y_ell, z_ell, alpha=0.2)  # alpha is the transparency

#plt.show()


#very inefficient way of doing this

