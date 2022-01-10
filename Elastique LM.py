import numpy as np
import matplotlib.pyplot as plt

from Elastique import pq2sig

E = 9 * 10**6 #*MPa
nu = 0.3
K = 0.01
G = 2 * 10**6 #*MPa
e0 = 0.80

pas = 100

H_hook = 1/E*np.array([ [1, -2*nu], 
                    [-nu, 1-nu]])

Me = np.array([ [1, 2],
                [2/3, -2/3]])

Mp = np.array([ [1/3, 2/3],
                [1, -1]])

def invPartiel(M):
    [[a,b], [c,d]] = M
    return np.array([[1/a, -a/b], [c/a, d-c*b/a]])

def H_Nlin(p):
    return np.array([[K/(1+e0)/p, 0], [0, 1/3*G]])

def indice_des_vides(ev, e0=e0):
    return -ev(1+e0)+e0

def Maitriser_la_contrainte():
    pas_sig = 10 #kN
    sigma_init = 500

    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([0,0])]
    PQ      = [np.array([0,0])]
    Epsilon_vp = [np.array([0,0])]

    while Sigma[-1][0]<sigma_init:

        
    return None

def Maitriser_contrainte_et_déplacement():
    pas_eps = 10 #sans unité
    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([0,0])]
    PQ      = [np.array([0,0])]
    Epsilon_vp = [np.array([0,0])]
    eps1=0
    while Sigma[-1][0] < 500 : 
        #calcule
        M = Me.inv @ H_hook @ Mp
        eps_1 = Epsilon[-1][0] + pas_eps
        sig_2 = 0
        sig_1, eps_2 = invPartiel(M)@np.array(eps_1, sig_2)

        #stock
        Sigma.append(np.array([sig_1, sig_2]))
        Epsilon.append(np.array([eps_1,eps_2]))
        PQ.append(Mp@Sigma[-1])
        Epsilon_vp.append(Me@Epsilon[-1])
    
    return Sigma, Epsilon, PQ, Epsilon_vp


def print(Sigma, Epsilon, PQ, Epsilon_vp):
    #TODO faire la fonction qui affiche tout bien
    return None

