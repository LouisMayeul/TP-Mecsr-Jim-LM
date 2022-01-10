import numpy as np
import matplotlib.pyplot as plt

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


def Maitriser_la_contrainte():
    #TODO faire la fonction while
    return None

def Maitriser_contrainte_et_déplacement():
    #TODO faire la fonction while avec la fonction invPartiel
    return None