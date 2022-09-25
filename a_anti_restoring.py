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
A = 10 #given in paper
m=23*(1.67e-27) #mass of 23Na
Lambda=589e-9 #wavelength of transition

N=2000  #number of atoms

z=np.arange(start=0,stop=200,step=0.1)  #so this is z, taking step = 0.1



# int a,b,c,d,e,f,g,h,i,j

# n_l = [[0]*len(z) for i in range(5)]     #initialized to 0, hope it doesnt affect the equations  #a 2D array of size len(z)x5
# n_u = [[0]*len(z) for i in range(5)]        #a 2D array of size len(z)x5
# np_array = np.array(n_l, dtype=np.int32)
# np_array = np.array(n_u, dtype=np.int32)

n_l = [1]*len(z)


# assuming each state of n_l (l=-2) is equally populated so each would contain N/len(z) atoms
# so
#print(len(z)) 

# for i in range(len(z)):
#     n_l[0][i] = 5   #all atoms are in ground state   so at each z value there are 5 atoms, so total atoms = 5*len(z) = 5*200 = 1000
# print(n_l)
# defining pairs that follow selection rule
M_u=[-2, -1, +0, +1, +2] # corresponding to F=2 (excited)
M_l=[-1, +0, +1]     # corresponding to F=1 (ground)

#     n_l[0][0]   n_l[0][1]  . . .   n_l[0][20mm]
     
#     n_l[1][0]   n_l[1][1]  . . .   n_l[1][20mm]    
#     n_l[2][0]   n_l[2][1]  . . .   n_l[2][20mm]   
#     n_l[3][0]   n_l[3][1]  . . .   n_l[1][20mm]  
#     n_l[4][0]   n_l[4][1]  . . .   n_l[1][20mm]

#      same with n_u

def g_l(l):
    return 0 #assumption
       
# Importing all the necessary libraries
def f( l,u,p ):
    return 1 # for now assuming all fractional strengths are 1

def detu( l,u,p ):
    return -tau # ( omega_p - omega( l,u) )  # The deturning

def zee_shift( l,u, g_u,z_val):
    A = 10   # given in paper
    num = g_u*M_u[u] - g_l(l)*M_l[l]
    return num*-2* A * z_val* mu_B/h_bar  # The Zeeman Shift

#We have also assumed that the branching ratio r(l,u) is 1 for all transitions
    

def R( l,u,p, z_val, g_u):
    u = l + 1 #for anti-restoring beam, deltaM = +1
    num = f( l,u,p )*s
    den = (1 + 4*( detu( l,u,p ) - 2*pi*0 - zee_shift(l,u, g_u,z_val) )**2)/(tau ** 2)
    #den = (1 + 4*( 1 - 2*pi*0 - z))/(tau ** 2)
    return tau*num/(2*den)
    

def r (l,u):
    return 1  #test





dipole_matrix_element  = 2.9883*e-29
d = dipole_matrix_element



#frac_s = [1/sqrt(30), 1/sqrt(10), 1/sqrt(5), 1/sqrt(3), 1/sqrt(2)]# fractional strength array
#for i in range(5):
 #   frac_s[i] = frac_s[i]*d
# ^ I think that these are the fractional strength values, just confirm it

# frac_s=[(3/4)*s, (1/2)*s, (1/4)*s, (3/4)*s, (1/2)*s, (1/4)*s]  #



# n_l= np.array(n_l)[indices.astype(int)]
# n_u= np.array(n_u)[indices.astype(int)]


#n_l[0][z], n_l[1][z], n_l[2][z], n_l[3][z], n_l[4][z], n_u[0][z], n_u[1][z], n_u[2][z], n_u[3][z], n_u[4][z] = symbols('n_l[0][z], n_l[1][z], n_l[2][z], n_l[3][z], n_l[4][z], n_u[0][z], n_u[1][z], n_u[2][z], n_u[3][z], n_u[4][z] ')

  
# defining equations
# eq1 = Eq( ( R(-2,-1,1,z)*(n_u[0][z] - n_l[0][z]) + tau*r(-2,-1)*n_u[0][z] ), 0)
# eq2 = Eq( ( R(-2,-1,1,z)*(n_u[1][z] - n_l[1][z]) + tau*r(-2,-1)*n_u[1][z] ), 0)
# eq3 = Eq( ( R(-2,-1,1,z)*(n_u[2][z] - n_l[2][z]) + tau*r(-2,-1)*n_u[2][z] ), 0)
# eq4 = Eq( ( R(-2,-1,1,z)*(n_u[3][z] - n_l[3][z]) + tau*r(-2,-1)*n_u[3][z] ), 0)
# eq5 = Eq( ( R(-2,-1,1,z)*(n_u[4][z] - n_l[4][z]) + tau*r(-2,-1)*n_u[4][z] ), 0)

# eq6 = Eq(  ( -( R(-2,-1,1,z)*(n_u[0][z] - n_l[0][z]) + tau*r(-2,-1)*n_u[0][z] ) -tau*n_u[0][z] ) , 0)
# eq7 = Eq(  ( -( R(-2,-1,1,z)*(n_u[1][z] - n_l[1][z]) + tau*r(-2,-1)*n_u[0][z] ) -tau*n_u[1][z] ) , 0)
# eq8 = Eq(  ( -( R(-2,-1,1,z)*(n_u[2][z] - n_l[2][z]) + tau*r(-2,-1)*n_u[0][z] ) -tau*n_u[2][z] ) , 0)
# eq9 = Eq(  ( -( R(-2,-1,1,z)*(n_u[3][z] - n_l[3][z]) + tau*r(-2,-1)*n_u[0][z] ) -tau*n_u[3][z] ) , 0)
# eq10 =Eq(  ( -( R(-2,-1,1,z)*(n_u[4][z] - n_l[4][z]) + tau*r(-2,-1)*n_u[0][z] ) -tau*n_u[4][z] ) , 0)

# eq11 = Eq( n_l[0] + n_l[1] + n_l[2] + n_l[3] + n_l[4] + n_u[0] + n_u[1] + n_u[2] + n_u[3] + n_u[4], N)   # sum of all n_l and n_u = N


# print(solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10 ), ( n_l[0], n_l[1], n_l[2], n_l[3], n_l[4], n_u[0], n_u[1], n_u[2], n_u[3], n_u[4] )))

#We must specify the initial conditions

# RHS = np.zeros( ( (len(z)*5) , 1) )  # column array of zeroes

# coeff_array = np.empty( ( (len(z) , 10 ) ) )

                       



#Trying to solve by Matrix method


#  A =  
#   -
#   |   C1,   C2,       0 ,     0 ,     0 ,     0 ,     0 ,      0 ,    0 ,      0      |
#      
#     n_l[0]  n_u[0]  n_l[1]  n_u[1]  n_l[2]  n_u[2]  n_l[3]  n_u[3]  n_l[4]  n_u[4]
#       
#       rows corresponding to different z positions follow
     
        
#   |  0,     0,      C3 ,     C4,      0,      0,      0,      0,      0,      0      |  
#   |  0,     0,       0,      0,       C5,     C6,      0,      0,      0,      0      |
#   |  0,     0,       0,      0,       0,      0,      C7,     C8,      0,      0      |
#   | 0,     0,       0,      0,       0,      0,      0,      0,      C9,     C10     |

# Now for the second set of equations obtained by putting  d(n_u)/dt = 0

#   | D1,     D2,       0,      0,       0,      0,      0,      0,      0,      0      |
#   | 0,       0,      D3,     D4,       0,      0,      0,      0,      0,      0      |
#   | 0,       0,       0,      0,       D5,     D6,      0,      0,      0,      0      | 
#   | 0,       0,       0,      0,       0,      0,      D7,     D8,      0,      0      |
#   | 0,       0,       0,      0,       0,      0,      0,      0,      D9,     D10     |







#al_z = [0 for i in range(2000)]  # empty acceleration array as a function of z


# for j in range(len(z)):
    
# al_z = -( R(-2,-1,1,z)*(n_u[0] - n_l[0]) + R(-2,-1,1,z)*(n_u[1] - n_l[1]) + R(-2,-1,1,z)*(n_u[2] - n_l[2]) + R(-2,-1,1,z)*(n_u[3] - n_l[3]) + R(-2,-1,1,z)*(n_u[4] - n_l[4]) )



def a_exp(g_u):
    a_z_val=[]
    for j in range(len(z)):
        a_z = (n_l[j]*R(0,2,1, (j*0.1) ,g_u)*tau) /( R(0,2,1, (j*0.1) ,g_u) + tau)   
        a_z_val.append(a_z)
    return a_z_val
     
al_z=[None]*11
for g_u in range(1,11):
    al_z[g_u]=a_exp(g_u/10)
 
       
# print(al_z[0])

# print (al_z) 


# population_arr = [ n_l[0],  n_u[0],  n_l[1],  n_u[1],  n_l[2],  n_u[2],  n_l[3],  n_u[3],  n_l[4],  n_u[4] ]  # array of population values for each z position

#   like  n_l[0][0]  n_l[0][1]  . . .  n_l[0][2000]  
#         n_u[0][0]  n_u[0][1]  . . .  n_u[0][2000] 
#         n_l[1][0]  n_l[1][1]  . . .  n_l[1][2000]
#         n_u[1][0]  n_u[1][1]  . . .  n_u[1][2000]
#         n_l[2][0]  n_l[2][1]  . . .  n_l[2][2000]
#         n_u[2][0]  n_u[2][1]  . . .  n_u[2][2000]
#         n_l[3][0]  n_l[3][1]  . . .  n_l[3][2000]
#         n_u[3][0]  n_u[3][1]  . . .  n_u[3][2000]
#         n_l[4][0]  n_l[4][1]  . . .  n_l[4][2000]
#         n_u[4][0]  n_u[4][1]  . . .  n_u[4][2000]

# print(al_z[3])         
  

plt.figure()
for g_u in range(1,11):
    plt.plot(z, al_z[g_u], label=f'g_u={g_u/10}', linewidth=2)
    
#plt.plot(z, al_z[i], label= f'g_u={i/10}', linewidth=1)


plt.title('gl=0; anti-restoring beam')
# plt.ylim([-600, 000])
# plt.xlim([0, 20])
plt.xlabel('z (mm)')
plt.ylabel('a_z (m/s^2')
plt.legend(loc='lower right', ncol=2, fontsize=7.5)
plt.show()

#fig = plt.figure(figsize=(2.79, 4.18))
#fig.show()  Dont know what this is for?







