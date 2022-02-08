import numpy as np
import matplotlib.pyplot as plt
import time

#from Elastique import pq2sig

E = 9 * 10**6 #*Pa
nu = 0.3
K = 0.01
G = 2 * 10**6 #*Pa
e0 = 0.80




Me = np.array([ [1, 2],
                [2/3, -2/3]])

Mp = np.array([ [1/3, 2/3],
                [1, -1]])

def invPartiel(M):
    [[a,b], [c,d]] = M
    return np.array([[1/a, -b/a], [c/a, d-c*b/a]])

def H(linéaire = True, p = 0):
    if linéaire:
        return 1/E*np.array([ [1, -2*nu], [-nu, 1-nu]])
    else:
        return np.linalg.inv(Me) @ np.array([[K/(1+e0)/p, 0], [0, 1/3/G]]) @ Mp

def indice_des_vides(ev, e0=e0):
    return -ev*(1+e0)+e0

def Maitriser_la_contrainte(pas_sig = 1000, linéaire = True):
    sigma_max = 500000
    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([10000,0])]
    PQ      = [Mp @ Sigma[-1]]
    Epsilon_vp = [np.array([0,0])]

    while Sigma[-1][0]<sigma_max:
        #calcul
        sig = np.array([Sigma[-1][0] + pas_sig, Sigma[-1][1]])
        #print(Sigma)
        eps = H(linéaire = linéaire,p = PQ[-1][0]) @ sig 
        print(eps)
        #Stock
        Sigma.append(sig)
        Epsilon.append(eps)
        PQ.append(Mp@sig)
        Epsilon_vp.append(Me@eps)

    return np.array(Sigma), np.array(Epsilon), np.array(PQ), np.array(Epsilon_vp)


def Maitriser_contrainte_et_déplacement(pas_eps = 1e-3, linéaire = True):
    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([10000,0])]
    PQ      = [Mp@Sigma[-1]]
    Epsilon_vp = [np.array([0,0])]
    i = 0
    while Sigma[-1][0] < 500000 : 
        print(Sigma[-1])
        #calcule
        M = H(linéaire = linéaire,p = PQ[-1][0])
        eps_1 = Epsilon[-1][0] + pas_eps
        print("eps1 = ", eps_1)
        sig_2 = 0
        sig_1, eps_2 = invPartiel(M) @ np.array([eps_1, sig_2])
        #stock
        Sigma.append(np.array([sig_1, sig_2]))
        Epsilon.append(np.array([eps_1,eps_2]))
        PQ.append(Mp@Sigma[-1])
        Epsilon_vp.append(Me@Epsilon[-1])
    
    return np.array(Sigma), np.array(Epsilon), np.array(PQ), np.array(Epsilon_vp)


def print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,title):
    #TODO faire la fonction qui affiche tout bien
    
    figure = plt.figure(figsize = (20,10)) # pour afficher les 4 courbes en même temps pour mieux comparer les différentes méthodes
    plt.figure(1)
    plt.suptitle(title)
    
    plt.subplot(2, 3, 1)
    plt.plot(Sigma[:,0],Sigma[:,1],color='#9F00EC')
    plt.title('$\sigma_1$ en fonction de $\sigma_3$')
    plt.ylabel('$\sigma_3$')
    plt.xlabel('$\sigma_1$')
    
    plt.subplot(2, 3, 2)
    plt.plot(PQ[:,0],PQ[:,1],color='green')
    plt.title("p en fonction de q")
    plt.ylabel('q')
    plt.xlabel('p')
    
    plt.subplot(2, 3, 3)
    plt.plot(PQ[10:,0],indice_des_vides(Epsilon_vp[10:,0]),color='blue')
    plt.title("p en fonction de e")
    plt.ylabel('e')
    plt.xlabel('p')
    
    plt.subplot(2, 3, 4)
    plt.plot(np.log(PQ[10:,0]),indice_des_vides(Epsilon_vp[10:,0]),color='red')
    plt.title("ln(p) en fonction de e")
    plt.ylabel('e')
    plt.xlabel('ln(p)')
    
    plt.subplot(2, 3, 5)
    plt.plot(Epsilon[10:,0],Epsilon_vp[10:,0],color='red')
    plt.title("$\epsilon_1$ en fonction de $\epsilon_v$")
    plt.ylabel('$\epsilon_3$')
    plt.xlabel('$\epsilon_1$')
    
    plt.subplot(2, 3, 6)
    plt.plot(Epsilon[10:,0],PQ[10:,1],color='green')
    plt.title("$\epsilon_1$ en fonction de q")
    plt.ylabel('q')
    plt.xlabel('$\epsilon_1$')
    plt.show()

    
    return None




 
Sigma, Epsilon, PQ, Epsilon_vp = Maitriser_la_contrainte()
print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,"Loi de Hooke - Contrainte")

Sigma, Epsilon, PQ, Epsilon_vp = Maitriser_contrainte_et_déplacement()
print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,"Loi de Hooke - Contrainte et Déplacement")

Sigma, Epsilon, PQ, Epsilon_vp = Maitriser_la_contrainte(linéaire = False)
print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,"Loi non Linéaire - Contrainte")

Sigma, Epsilon, PQ, Epsilon_vp = Maitriser_contrainte_et_déplacement(linéaire = False)
print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,"Loi non Linéaire - Contrainte et Déplacement")
