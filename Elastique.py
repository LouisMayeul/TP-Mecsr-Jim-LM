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

E = 9 * 10**6 #*MPa
nu = 0.3
K = 0.01
G = 2 * 10**6 #*MPa
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

def sig2pq(sigma):
    """ transforme le vecteur contenant sigma1 et sigma3 
    en u vecteur contenant p et q"""
    return np.array([(sigma[0] + 2*sigma[1])/3, sigma[0] - sigma[1]])

def pq2sig(p,q):
    """ transforme p et q
    en u vecteur contenant sigma1 et sigma3"""
    return np.array([p + 2*q/3, p - q/3 ])




def epsilon2eveq(epsilon):
    """ donne le vecteur de epsilon_v et epsilon_q
    en fonction du vecteur epsilon """
    return np.array([epsilon[0] + 2*epsilon[1], 2*(epsilon[0] - epsilon[1])/3])

def eveq2epsilon(evq):
    """ donne le vecteur epsilon
    en fonction du vecteur epsilon_v et epsilon_q """
    return np.array([evq[0]/3 + evq[1], evq[0]/3 - evq[1]/2])


############################################################
######################### Lois : ###########################
############################################################


def loi_non_lineraire_de2pq(dev,deq, p_avt):
    """ permet d'optenir p et q grâce à la loi non linéaire donnée
    en fonction de d(epsilon_v) et d(epsilon_q)"""
    return np.array([(1+e0)/K*p_avt*dev, 3*G*deq])

def loi_Hook_de2pq(dev,deq):
    """ permet d'optenir p et q grâce à la loi de Hook
    en fonction de d(epsilon_v) et d(epsilon_q)"""
    #TODO à faire
    return np.array([0,0])


def loi_non_lineraire_pq2de(dp,dq, p_avt):
    """ permet d'optenir d(epsilon_v) et d(epsilon_q) 
    grâce à la loi non linéaire donnée
    en fonction de p et q """
    return np.array([(K*dp)/((1+e0)*p_avt), dq/(3*G)])

def loi_Hook_pq2de(dev,deq):
    """ permet d'optenir d(epsilon_v) et d(epsilon_q) 
    grâce à la loi de Hook
    en fonction de p et q """
    #TODO à faire
    return np.array([0,0])

############################################################
####################### Programme : ########################
############################################################


def calcul(sigma1,sigma3,epsilon):
    
    while sigma1[-1] < sigma_1f :
        
        sigma_n = np.array([sigma1[-1],sigma3[-1]]) #on récupère la valeur actuelle de sigma
        sigma1.append(sigma1[-1] + pas) #on augmente sigma1 pas à pas
        sigma3.append(sigma3[-1]) #on ne change pas la valeur de sigma3
        sigma_n1 = np.array([sigma1[-1],sigma3[-1]]) #on récupère le sigma à l'intant infinitésimal suivant
        dsigma = sigma_n1 - sigma_n #on calcul la variation de sigma
        
        p_n1 = sig2pq(sigma_n1)[0] #on calcule p à l'intant suivant
        d = sig2pq(dsigma)
        dp = d[0] #on calcule la variation de p dp en fonction de la variation de sigma dsigma
        dq = d[1] #on calcule la variation de q dq en fonction de la variation de sigma dsigma
        
        devq = loi_non_lineraire_pq2de(dp,dq,p_n1) #de même pour la variation de epsilon_v
        
        #on calcule les variations de epsilon_1 et epsilon_3 en fonction des epsilon_v et epsilon_q
        d_epsilon = epsilon2eveq(devq)
        
        #on rajoute les valeurs de epsilon dans la suite
        epsilon.append([epsilon[-2][0]+d_epsilon[0],epsilon[-2][1]+d_epsilon[1]])
        
    return sigma1,sigma3,epsilon

############################################################
####################### Question 1 : #######################
############################################################

sigma1 = [sigma_1i]
sigma3 = [sigma_3i]
epsilon = [ [0,0],
            [0,0]]

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
    P.append(sig2pq([s1[k],s3[k]])[0])
    Q.append(sig2pq([s1[k],s3[k]])[1])
    Eps.append(e[k][0] + 2*e[k][1])  #! Ca c'est ev
    E1.append(e[k][0])
    Ev.append(epsilon2eveq(e[k]))

### Affichage

plt.plot(s1,s3,label='$\sigma_1$ en fonction de $\sigma_3$')  #! sigma1 est nul
plt.legend()
plt.grid(True)
plt.show()

plt.plot(P,Q,label='p en fonction de q') #! linéaire c'est 1/3 c'est logique
plt.legend()
plt.grid(True)
plt.show()

plt.plot(P,Eps,label='p en fonction de e') #! PAS E VRAIMENT...pk pas
plt.legend()
plt.grid(True)
plt.show()

plt.plot(np.log(P),Eps,label='ln(p) en fonction de e') #linéaire mais ok
plt.legend()
plt.grid(True)
plt.show()

plt.plot(E1,Ev,label='$\epsilon_1$ en fonction de $\epsilon_v$') #pk pas
plt.legend()
plt.grid(True)
plt.show()

plt.plot(E1,Q,label='$\epsilon_1$ en fonction de q') #pk pas
plt.legend()
plt.grid(True)
plt.show()

############################################################
####################### Question 2 : #######################
############################################################


def Calcul2(sigma_n,sigma_f,defo_n):
    """ fait le calcul avec epsilon fixé"""
    print("calcul 2")
    print(defo_n)
    #On différencie les epsilones v et q
    dev, deq =defo_n[-1]-defo_n[-2]

    #ce qui nous donne les variations de p et qavec la loi non linéaire 
    dp = (1+e0)/K*sig2pq(sigma_n[-1])[0]*dev
    dq = 3*G*deq

    #On peut récuperer les valeurs de p et q    
    p,q = sig2pq(sigma_n[-1])+ np.array([dp,dq])
    #qui sont transformé en sigma
    sigma_n.append(pq2sig(p,q))

    while sigma_n[-1][0] < sigma_f[0]:
        #On fait cela tant que sigma1 est inf à la valeur seuil
        dev, deq =defo_n[-1]-defo_n[-2]
        defo_n .append(defo_n[-1] - np.array([pas, 0]))

        dp = (1+e0)/K*sig2pq(sigma_n[-1])[0]*dev
        dq = 3*G*deq

        p,q = sig2pq(sigma_n[-1])+ np.array([dp,dq])

        sigma_n.append(sigma_n[-1] + pq2sig(p,q))
        print(sigma_n[-1])
    return sigma_n, defo_n


def Q2():
    print("q 2")
    sigma_n = [[10,0]]
    sigma_f = [500e3, 0]

    defo_n = [np.array([0,0]),
              np.array([10,0])]
    defo_f = []

    sigma_n, defo_n=Calcul2(sigma_n,sigma_f, defo_n)
    print(sigma_n, defo_n)

#Q2()

    




















