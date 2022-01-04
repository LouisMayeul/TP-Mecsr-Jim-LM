# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 23:35:03 2022

@author: Jim
"""
import numpy as np
import matplotlib.pyplot as plt

############################################################
######################## Données : #########################
############################################################

E = 9 * 10**6
nu = 0.3
K = 0.01
G = 2 * 10**6
e0 = 0.80

pas = 100

### Question 1 
sigma_1i = 10 * 10**3
sigma_2i = 0
sigma_3i = 0

sigma_1f = 500 * 10**3
sigma_2f = 0
sigma_3f = 0



############################################################
####################### Fonctions : ########################
############################################################

def p(sigma):
    return (sigma[0] + 2*sigma[1])/3

def q(sigma):
    return sigma[0] - sigma[1]

def sigma1(p,q):
    return p + 2*q/3

def sigma3(p,q):
    return p - q/3



def ev(epsilon):
    return epsilon[0] + 2*epsilon[1]

def eq(epsilon):
    return 2*(epsilon[0] - epsilon[1])/3

def e1(ev,eq):
    return ev/3 + eq

def e3(ev,eq):
    return ev/3 - eq/2

############################################################
####################### Programme : ########################
############################################################

sigma1 = [sigma_1i]
sigma3 = [sigma_3i]
epsilon = [[0,0],[0,0]]


def calcul(sigma1,sigma3,epsilon):
    
    while sigma1[-1] < sigma_1f :
        
        sigma_n = [sigma1[-1],sigma3[-1]]
        sigma1.append(sigma1[-1] + pas)
        sigma3.append(sigma3[-1])
        sigma_n1 = [sigma1[-1],sigma3[-1]]
        dsigma = [sigma_n1[0] - sigma_n[0],sigma_n1[1] - sigma_n[1]]
        
        p_n = p(sigma_n1)
        dp = p(dsigma)
        dq = q(dsigma)
        
        dev = (K*dp)/((1+e0)*p_n)
        deq = dq/(3*G)
        
        de1 = e1(dev,deq)
        de3 = e3(dev,deq)
        
        epsilon.append([epsilon[-2][0]+de1,epsilon[-2][1]+de3])
        
    return sigma1,sigma3,epsilon


s1,s3,e = calcul(sigma1, sigma3, epsilon)
s1.pop(-1)
s3.pop(-1)
e.pop(-1)
P = []
Q = []
Eps = []
E1 = []
Ev = []


for k in range(len(s1)):
    P.append(p([s1[k],s3[k]]))
    Q.append(q([s1[k],s3[k]]))
    Eps.append(e[k][0] + 2*e[k][1])
    E1.append(e[k][0])
    Ev.append(ev(e[k]))

print(Ev)



plt.plot(s1,s3,label='$\sigma_1$ en fonction de $\sigma_3$')  #Ne sigma1 est nul
plt.legend()
plt.grid(True)
plt.show()
plt.plot(P,Q,label='p en fonction de q')
plt.legend()
plt.grid(True)
plt.show()
plt.plot(P,Eps,label='p en fonction de e')
plt.legend()
plt.grid(True)
plt.show()
plt.plot(np.log(P),Eps,label='ln(p) en fonction de e')
plt.legend()
plt.grid(True)
plt.show()
plt.plot(E1,Ev,label='$\epsilon_1$ en fonction de $\epsilon_v$')
plt.legend()
plt.grid(True)
plt.show()
plt.plot(E1,Q,label='$\epsilon_1$ en fonction de q')

plt.legend()
plt.grid(True)
plt.show()





















