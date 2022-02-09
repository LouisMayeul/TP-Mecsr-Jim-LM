from Elastique_LM import *
import numpy as np

E = 9 * 1e6 #*Pa
nu = 0.3
K = 0.01
G = 2 * 10**6 #*Pa
e0 = 0.80

Me = np.array([ [1, 2],
                [2/3, -2/3]])

Mp = np.array([ [1/3, 2/3],
                [1, -1]])


M = 1,0
p0_init = 100*1e3
lamdba = 0.2


pv = 1e-3 #Petite valeur pour égalité à 0



def H_tilte_plastique(pq, p0):
    p,q=pq
    a = M**2*(2*p-p0)/(M**2*p*p0*(1+e0)/(lamdba -K)) + K(1+e0)*p
    b = 2*q/(M**2*p*p0*(1+e0)/(lamdba-K))
    c = 2*q*M**2*(2*p-p0)/(M**4*p*p0*(1+e0)/(lamdba-K)*(2*p-p0))
    d = 4*q**2/(M**4*p*p0*(1+e0)/(lamdba-K)*(2*p-p0)) + 1/3/G
    return np.array([[a,b],[c,d]])

def f(pq,p0):
    """définit la surface de charge"""
    p,q=pq
    return q**2-M**2*p*(p0-p)

def display_surface_de_charge(p0, nb_pt = 100):
    """définit la surface de charge qcharge en fonction de p"""
    PQ=[]
    for i in range(nb_pt):
        p = i*p0/nb_pt
        q = np.sqrt(M**2*p*(p0-p))
        PQ.append(np.array([p,q]))
    return np.array(PQ)



def essai_oedométrique(pas_eps = 1e-4):

    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([10* 1e3,5* 1e3])]
    PQ      = [Mp @ Sigma[-1]]
    Epsilon_vp = [Me@Epsilon[-1]]
    p0 = p0_init

    d_eps = np.array([pas_eps,0])
    if f(PQ[-1], p0) <= pv:
        H = np.linalg.inv(Me) @ H(linéaire = True, PQ[-1][0]) @ Mp
        d_sig = np.linalg.inv(H)@
    else:
        H = np.linalg.inv(Me) @ H_tilte_plastique(PQ[-1],p0) @ Mp
        pass


