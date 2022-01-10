import numpy as np
import matplotlib.pyplot as plt

from Elastique import pq2sig

E = 9 * 10**6 #*MPa
nu = 0.3
K = 0.01
G = 2 * 10**6 #*MPa
e0 = 0.80

pas = 100

Me = np.array([ [1, 2],
                [2/3, -2/3]])

Mp = np.array([ [1/3, 2/3],
                [1, -1]])

def invPartiel(M):
    [[a,b], [c,d]] = M
    return np.array([[1/a, -a/b], [c/a, d-c*b/a]])

def H(linéaire = True, p = 0):
    if linéaire:
        return 1/E*np.array([ [1, -2*nu], [-nu, 1-nu]])
    else:
        return np.array([[K/(1+e0)/p, 0], [0, 1/3*G]])

def indice_des_vides(ev, e0=e0):
    return -ev(1+e0)+e0

def Maitriser_la_contrainte(linéaire = True):
    pas_sig = np.array([10,0]) #kN
    sigma_init = 500

    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([0,0])]
    PQ      = [np.array([0,0])]
    Epsilon_vp = [np.array([0,0])]

    while Sigma[-1][0]<sigma_init:
        #calcul
        sig = Sigma[-1] + pas_sig
        eps = H(linéaire, PQ[-1][0])@sig 
    
        #Stock
        Sigma.append(sig)
        Epsilon.append(eps)
        PQ.append(Mp@sig)
        Epsilon_vp.append(Me@eps)

    return Sigma, Epsilon, PQ, Epsilon_vp


#pregunta para théo : el metro de los muetos, ?quitaron los asientos?
def Maitriser_contrainte_et_déplacement(linéaire = True):
    pas_eps = 10 #sans unité
    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([0,0])]
    PQ      = [np.array([0,0])]
    Epsilon_vp = [np.array([0,0])]
    eps1=0
    while Sigma[-1][0] < 500 : 
        #calcule
        M = Me.inv @ H(linéaire, PQ[-1][0]) @ Mp
        eps_1 = Epsilon[-1][0] + pas_eps
        sig_2 = 0
        sig_1, eps_2 = invPartiel(M)@np.array(eps_1, sig_2)
        print(sig_1)
        #stock
        Sigma.append(np.array([sig_1, sig_2]))
        Epsilon.append(np.array([eps_1,eps_2]))
        PQ.append(Mp@Sigma[-1])
        Epsilon_vp.append(Me@Epsilon[-1])
    
    return Sigma, Epsilon, PQ, Epsilon_vp


def print_graphe(Sigma, Epsilon, PQ, Epsilon_vp):
    #TODO faire la fonction qui affiche tout bien
    return None


print(Maitriser_contrainte_et_déplacement())
