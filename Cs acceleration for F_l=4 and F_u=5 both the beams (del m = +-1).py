from cmath import sqrt
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# Importing Lib_s
# Defining Parameters

Gamma=2*(22/7)*(5.33e6)
Delta=-1*Gamma
h=6.626e-34
m=133*(1.67e-27)
Lambda=852e-9
v=0
Mu_B=9.274e-24
A=0.1
H_bar=1.055e-34
Pi=22/7


# Playing with numbers to reduce complexicity. Refer to final formula

p=(-1*h*Gamma)/(2*m*Lambda)
# p is propotional to force


q=((2*Mu_B*A)/(H_bar*Gamma))*(0.001)
# q is the factor in zeeman shift. 0.001 is just changing scale.

#these are the constant factors

g_l=-1
# g_l is the g factor for lower level

# defining pairs that follow selection rule
M_u=[-5, -4, -3, -2, -1, +0, +1, +2, +3, +5, +4, +3, +2, +1, +0, -1, -2, -3]
M_l=[-4, -3, -2, -1, +0, +1, +2, +3, +4, +4, +3, +2, +1, +0, -1, -2, -3, -4]


# assume f to be 5 for all
f=[1/sqrt(2), 1/sqrt(5/2), 1/sqrt(45/14), 1/sqrt(30/7), 1/sqrt(6), 1/sqrt(9), 1/sqrt(15), 1/sqrt(30), 1/sqrt(90), 1/sqrt(2), 1/sqrt(5/2), 1/sqrt(45/14), 1/sqrt(30/7), 1/sqrt(6), 1/sqrt(9), 1/sqrt(15), 1/sqrt(30), 1/sqrt(90)]

# defining z axis
z=np.arange(start=0,stop=20,step=0.01)

# func1 is for laser beam going to Right/Left(Have to check methamatically. Doesnt matter for now if we are getting correct behavior)
def func1(k,i):
    return (f[k]/(1+f[k]+4*((1+q*z*((i/10)*M_u[k]-g_l*M_l[k]))**2)))
        
        
# func2 is for laser beam going to Right/Left(Have to check methamatically. Doesnt matter for now if we are getting correct behavior)          
def func2(k,i):
    return (f[k+9]/(1+f[k+9]+4*((1+q*z*((i/10)*M_u[k+9]-g_l*M_l[k+9]))**2)))


# func3 is force due to both
def func3(k, i):
    return (func1(k, i)-func2(k, i))*(p/10000)

# Defining force for all nine transition

def a_exp(i):
    a_z_val=[]
    for j in range(len(z)):
        a_z=(func3(0,i)+func3(1,i)+func3(2,i)+func3(3,i)+func3(4,i)+func3(5,i)+func3(6,i)+func3(7,i)+func3(8,i))
        a_z_val.append(a_z)
    return a_z

# Calculating and combining all graphs for each g_u value(Colorful graphs)

al_z=[None]*11
for i in range(1,11):
    al_z[i]=a_exp(i)
    
plt.figure()
for i in range(1,11):
    plt.plot(z, al_z[i], label=f'g_u={i/10}', linewidth=2)
plt.title('gl=-1')
# plt.ylim([-4, 0.5])
# plt.xlim([0, 10])
plt.xticks(fontsize=15, rotation=90)
plt.yticks(fontsize=15, rotation=0)
plt.xlabel('z (mm)', fontsize=15)
plt.ylabel('$ a_z (10^4 m/s^2) $', fontsize=15)
plt.legend(loc='lower right', ncol=2, fontsize=7.5)
fig = plt.figure(figsize=(2.79, 4.18))
plt.show()


