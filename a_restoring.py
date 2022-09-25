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
#from tkinter import N
from sympy import symbols, Eq, solve     # for solving equations



mu_B = 9.2740100783e-24
h_bar = 1.055e-34 # reduced Planck's constant
pi = 22.0/7.0     # pi
tau = 61.542e6     # Natural Line Width (F W H M)                                           
s = 1     #  given in paper
A = 1000 #given in paper
m=23*(1.67e-27) #mass of 23Na
Lambda=589e-9 #wavelength of transition


N=1000 #number of atoms

z=np.arange(start=0,stop=100,step=0.1)  #so this is z, taking step = 0.1


n_l = [1]*len(z)


# assuming each state of n_l (l=-2) is equally populated so each would contain N/len(z) atoms

# defining pairs that follow  del m  = -1 selection rule
M_u=[-2, -1, +0] # corresponding to F=2 (excited)
M_l=[-1, +0, +1]     # corresponding to F=1 (ground)



def g_l(l):
    return 0#assumption
       
# Importing all the necessary libraries
def f( l,u,p ):
    return 1 # for now assuming all fractional strengths are 1

def detu( l,u,p ):
    return -tau # ( omega_p - omega( l,u) )  # The deturning

def zee_shift( l,u, g_u,z_val):
    A = 1000   # given in paper
    num = g_u*M_u[u] - g_l(l)*M_l[l]
    return num*-2* A * z_val* mu_B/h_bar  # The Zeeman Shift

#We have also assumed that the branching ratio r(l,u) is 1 for all transitions
    

def R( l,u,p, z_val, g_u):
    # u = l - 1 #for restoring beam, deltaM = -1
    num = f( l,u,p )*s
    den = (1 + 4*( detu( l,u,p ) - 2*pi*0 - zee_shift(l,u, g_u,z_val) )**2/(tau ** 2) )
    #den = (1 + 4*( 1 - 2*pi*0 - z))/(tau ** 2)
    return tau*num/(2*den)
    

def r (l,u):
    return 1  #test





dipole_matrix_element  = 2.9883*e-29
d = dipole_matrix_element



frac_s = [1/sqrt(30), 1/sqrt(10), 1/sqrt(5), 1/sqrt(3), 1/sqrt(2)]# fractional strength array
for i in range(5):
    frac_s[i] = frac_s[i]*d
# ^ I think that these are the fractional strength values, just confirm it


# for j in range(len(z)):
    
# al_z = -( R(-2,-1,1,z)*(n_u[0] - n_l[0]) + R(-2,-1,1,z)*(n_u[1] - n_l[1]) + R(-2,-1,1,z)*(n_u[2] - n_l[2]) + R(-2,-1,1,z)*(n_u[3] - n_l[3]) + R(-2,-1,1,z)*(n_u[4] - n_l[4]) )



def a_exp(g_u):
    a_z_val=[]
    for j in range(len(z)):
        a_z = -(n_l[j]*R(0,0,-1, (j*10e-3) ,g_u)*tau) /( R(0,0,-1, (j*10e-3) ,g_u) + tau)   
        a_z_val.append(a_z)
    return a_z_val
     
al_z=[None]*11
for g_u in range(1,11):
    al_z[g_u]=a_exp(g_u/10)
 
       


plt.figure()
for g_u in range(1,11):
    plt.plot(z, al_z[g_u], label=f'g_u={g_u/10}', linewidth=2)
    
#plt.plot(z, al_z[i], label= f'g_u={i/10}', linewidth=1)


plt.title('gl=0; restoring')
# plt.ylim([-600, 000])
# plt.xlim([0, 20])
plt.xlabel('z (mm)')
plt.ylabel('a_z (m/s^2')
plt.legend(loc='lower right', ncol=2, fontsize=7.5)
plt.show()

#fig = plt.figure(figsize=(2.79, 4.18))
#fig.show()  Dont know what this is for?







