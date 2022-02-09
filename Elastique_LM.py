import numpy as np
import matplotlib.pyplot as plt
import time

#from Elastique import pq2sig

E = 9 * 1e6 #*Pa
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

def H_elastique(linéaire = True, p = 0):
    if linéaire:
        return 1/E*np.array([ [1, -2*nu], [-nu, 1-nu]])
    else:
        return np.linalg.inv(Me) @ np.array([[K/(1+e0)/p, 0], [0, 1/3/G]]) @ Mp

def indice_des_vides(ev, e0=e0):
    return -ev*(1+e0)+e0

def Maitriser_la_contrainte(pas_sig = 1* 1e3, linéaire = True):
    sigma_max = 500 * 1e3
    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([10* 1e3,0])]
    PQ      = [Mp @ Sigma[-1]]
    Epsilon_vp = [Me@Epsilon[-1]]

    d_sig = np.array([pas_sig,0])

    while Sigma[-1][0]<sigma_max:
        #calcul

        #print(Sigma)
        d_eps = H_elastique(linéaire = linéaire,p = PQ[-1][0]) @ d_sig 
    
        #Stock
        Sigma.append(Sigma[-1] + d_sig)
        Epsilon.append(Epsilon[-1] + d_eps)
        PQ.append(Mp@Sigma[-1])
        Epsilon_vp.append(Me@Epsilon[-1])

    return np.array(Sigma), np.array(Epsilon), np.array(PQ), np.array(Epsilon_vp)


def Maitriser_contrainte_et_déplacement(pas_eps = 1e-4, linéaire = True):
    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([10* 1e3,0])]
    PQ      = [Mp@Sigma[-1]]
    Epsilon_vp = [np.array([0,0])]

    d_eps_1 = pas_eps
    d_sig_2 = 0

    while Sigma[-1][0] < 500 * 1e3: 
        #calcule
        d_sig_1, d_eps_2 = invPartiel(H_elastique(linéaire = linéaire,p = PQ[-1][0])) @ np.array([d_eps_1, d_sig_2])

        #stock
        Sigma.append(Sigma[-1]+np.array([d_sig_1, d_sig_2]))
        Epsilon.append(Epsilon[-1]+ np.array([d_eps_1,d_eps_2]))
        PQ.append(Mp@Sigma[-1])
        Epsilon_vp.append(Me@Epsilon[-1])
    
    return np.array(Sigma), np.array(Epsilon), np.array(PQ), np.array(Epsilon_vp)


def print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,title):
    
    figure = plt.figure(figsize = (20,10)) # pour afficher les 4 courbes en même temps pour mieux comparer les différentes méthodes
    plt.figure(1)
    plt.suptitle(title)
    plt.subplots_adjust(wspace=0.3, hspace=.3)
    
    plt.subplot(2, 3, 1)
    plt.plot(Sigma[1:,0],PQ[1:,1],color='#9F00EC')
    plt.title(' $q$ en fonction de $\sigma_1$')
    plt.ylabel('$q$')
    plt.xlabel('$\sigma_1$')
    
    plt.subplot(2, 3, 2)
    plt.plot(PQ[1:,0],PQ[1:,1],color='green')
    plt.title("q en fonction de p")
    plt.ylabel('q')
    plt.xlabel('p')
    
    plt.subplot(2, 3, 3)
    plt.plot(PQ[1:,0],indice_des_vides(Epsilon_vp[1:,0]),color='blue')
    plt.title("e en fonction de p")
    plt.ylabel('e')
    plt.xlabel('p')
    
    plt.subplot(2, 3, 4)
    plt.plot(np.log(PQ[1:,0]),indice_des_vides(Epsilon_vp[1:,0]),color='red')
    plt.title("e en fonction de ln(p) ")
    plt.ylabel('e')
    plt.xlabel('ln(p)')
    
    plt.subplot(2, 3, 5)
    plt.plot(Epsilon[1:,0],Epsilon_vp[1:,0],color='red')
    plt.title("$\epsilon_v$ en fonction de $\epsilon_1$")
    plt.ylabel('$\epsilon_3$')
    plt.xlabel('$\epsilon_1$')
    
    plt.subplot(2, 3, 6)
    plt.plot(Epsilon[1:,0],PQ[1:,1],color='green')
    plt.title(" q en fonction de $\epsilon_1$")
    plt.ylabel('q')
    plt.xlabel('$\epsilon_1$')
    
    #plt.savefig(title+'.png')
    plt.show()

    
    return None


def partie1():
    taille=[]
    
    Sigma, Epsilon, PQ, Epsilon_vp = Maitriser_la_contrainte()
    print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,"Loi de Hooke - Contrainte")
    taille.append(len(Sigma))

    Sigma, Epsilon, PQ, Epsilon_vp = Maitriser_contrainte_et_déplacement(linéaire = True)
    print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,"Loi de Hooke - Contrainte et Déplacement")
    taille.append(len(Sigma))
    Sigma, Epsilon, PQ, Epsilon_vp = Maitriser_la_contrainte(linéaire = False)
    print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,"Loi non Linéaire - Contrainte")
    taille.append(len(Sigma))
    Sigma, Epsilon, PQ, Epsilon_vp = Maitriser_contrainte_et_déplacement(linéaire = False)
    print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,"Loi non Linéaire - Contrainte et Déplacement")
    taille.append(len(Sigma))

    print(taille)



#####Partie 2#############################
M = 1.0
p0_init = 100*1e3
lamdba = 0.2

pv = 1e2 #Petite valeur pour égalité à 0

def H_plastique(pq, p0):

    p,q=pq
    a = M**2*(2*p-p0)/(M**2*p*p0*(1+e0)/(lamdba -K)) + K/(1+e0)/p
    b = 2*q/(M**2*p*p0*(1+e0)/(lamdba-K))
    c = 2*q*M**2*(2*p-p0)/(M**4*p*p0*(1+e0)/(lamdba-K)*(2*p-p0))
    d = 4*q**2/(M**4*p*p0*(1+e0)/(lamdba-K)*(2*p-p0)) + 1/3/G

    return np.linalg.inv(Me) @ np.array([[a,b],[c,d]]) @ Mp

def Mat(pq, pr0):
    pr, q = pq
    Mat =  np.array([[M ** 2 * (2 * pr - pr0) / pr / pr0 / M ** 2 / (1 + e0) * (lamdba - K) + K / pr / (1 + e0),
                      2 * q / M ** 2 / pr / pr0 / (1 + e0) * (lamdba - K)],
                     [2 * q * M ** 2 * (2 * pr - pr0) / M ** 4 / pr / pr0 / (1 + e0) / (2 * pr - pr0) * (lamdba - K),
                      4 * q ** 2 / M ** 4 / pr / pr0 / (1 + e0) / (2 * pr - pr0) * (lamdba - K)+1/3/G]])

    return np.linalg.inv(Me) @ Mat @ Mp


def f(pq,p0):
    """définit la surface de charge"""
    p,q=pq
    return q**2-M**2*p*(p0-p)

def display_surface_de_charge(p0, nb_pt = 10000):
    """définit la surface de charge qcharge en fonction de p"""
    PQ=[]
    for i in range(nb_pt):
        p = i*p0/nb_pt
        q = np.sqrt(M**2*p*(p0-p))
        PQ.append(np.array([p,q]))
    return np.array(PQ)

def p0_uptudate(pq):
    p,q = pq
    return q**2/M**2/p + p

def essai_oedométrique(pas_eps = 1e-4):

    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([10* 1e3,5* 1e3])]
    PQ      = [Mp @ Sigma[-1]]
    Epsilon_vp = [Me@Epsilon[-1]]
    p0 = p0_init
    F = [f(PQ[-1], p0)]

    d_eps = np.array([pas_eps,0])


    while Epsilon[-1][0] < .2:

        if F[-1] <= -pv:
            d_sig = np.linalg.inv(H_elastique(linéaire = False, p = PQ[-1][0])) @ d_eps

            Sigma.append(Sigma[-1] + d_sig)
            Epsilon.append(Epsilon[-1] + d_eps)
            PQ.append(Mp@Sigma[-1])
            Epsilon_vp.append(Me@Epsilon[-1])
            F.append(f(PQ[-1], p0))

        else:

            d_sig = np.linalg.inv(H_plastique(PQ[-1],p0)) @ d_eps

            Sigma.append(Sigma[-1] + d_sig)
            Epsilon.append(Epsilon[-1] + d_eps)
            PQ.append(Mp@Sigma[-1])
            Epsilon_vp.append(Me@Epsilon[-1])

            p0 = p0_uptudate(PQ[-1])
            F.append(f(PQ[-1], p0))
    return np.array(Sigma), np.array(Epsilon), np.array(PQ), np.array(Epsilon_vp), F,p0



def essai_isotrope_draine(pas_sig=1e3):
    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([10* 1e3,10* 1e3])]
    PQ      = [Mp @ Sigma[-1]]
    Epsilon_vp = [Me@Epsilon[-1]]
    p0 = p0_init
    F = [f(PQ[-1], p0)]

    d_sig = np.array([pas_sig,pas_sig])


    while Sigma[-1][0] < 250 * 1e3:

        if F[-1] <= -pv:
            d_eps = H_elastique(linéaire = False, p = PQ[-1][0]) @ d_sig

            Sigma.append(Sigma[-1] + d_sig)
            Epsilon.append(Epsilon[-1] + d_eps)
            PQ.append(Mp@Sigma[-1])
            Epsilon_vp.append(Me@Epsilon[-1])
            F.append(f(PQ[-1], p0))

        else:

            d_eps = H_plastique(PQ[-1],p0) @ d_sig

            Sigma.append(Sigma[-1] + d_sig)
            Epsilon.append(Epsilon[-1] + d_eps)
            PQ.append(Mp@Sigma[-1])
            Epsilon_vp.append(Me@Epsilon[-1])

            p0 = p0_uptudate(PQ[-1])
            F.append(f(PQ[-1], p0))
    return np.array(Sigma), np.array(Epsilon), np.array(PQ), np.array(Epsilon_vp), F, p0


def essai_triaxial_draine(pas_sig=1e3):
    Sigma, Epsilon, PQ, Epsilon_vp, F, p0 = essai_isotrope_draine()

    Sigma = Sigma.tolist()
    Epsilon = Epsilon.tolist()
    PQ=PQ.tolist()
    Epsilon_vp=Epsilon_vp.tolist()

    d_sig = np.array([-pas_sig,-pas_sig])
    T = [len(Sigma)]

    #détente isotropique
    while Sigma[-1][0] > 200 * 1e3:

        d_eps = H_elastique(linéaire = False, p = PQ[-1][0]) @ d_sig

        Sigma.append(Sigma[-1] + d_sig)
        Epsilon.append(Epsilon[-1] + d_eps)
        PQ.append(Mp@Sigma[-1])
        Epsilon_vp.append(Me@Epsilon[-1])
        F.append(f(PQ[-1], p0))
    T.append(len(Sigma))

    d_sig = np.array([pas_sig,0])
    while PQ[-1][1] < M*PQ[-1][0] - pv:

        if F[-1] <= -pv:
            d_eps = H_elastique(linéaire = False, p = PQ[-1][0]) @ d_sig

            Sigma.append(Sigma[-1] + d_sig)
            Epsilon.append(Epsilon[-1] + d_eps)
            PQ.append(Mp@Sigma[-1])
            Epsilon_vp.append(Me@Epsilon[-1])
            F.append(f(PQ[-1], p0))

        else:

            d_eps = H_plastique(PQ[-1],p0) @ d_sig

            Sigma.append(Sigma[-1] + d_sig)
            Epsilon.append(Epsilon[-1] + d_eps)
            PQ.append(Mp@Sigma[-1])
            Epsilon_vp.append(Me@Epsilon[-1])

            p0 = p0_uptudate(PQ[-1])
            F.append(f(PQ[-1], p0))
    T.append(len(Sigma))
    return np.array(Sigma), np.array(Epsilon), np.array(PQ), np.array(Epsilon_vp), F,p0



def essai_triaxial_draine_surconso(pas_sig=1e3, pas_esp=1e-3):
    Sigma, Epsilon, PQ, Epsilon_vp, F, p0 = essai_isotrope_draine()

    Sigma = Sigma.tolist()
    Epsilon = Epsilon.tolist()
    PQ=PQ.tolist()
    Epsilon_vp=Epsilon_vp.tolist()

    d_sig = np.array([-pas_sig,-pas_sig])
    T = [len(Sigma)]

    #détente isotropique
    while Sigma[-1][0] > 50 * 1e3:

        d_eps = H_elastique(linéaire = False, p = PQ[-1][0]) @ d_sig

        Sigma.append(Sigma[-1] + d_sig)
        Epsilon.append(Epsilon[-1] + d_eps)
        PQ.append(Mp@Sigma[-1])
        Epsilon_vp.append(Me@Epsilon[-1])
        F.append(f(PQ[-1], p0))
    T.append(len(Sigma))

    d_eps_1 = pas_esp
    d_sig_3 = 0

    while PQ[-1][1] < M*PQ[-1][0] - pv:

        if F[-1] <= -pv:

            d_sig_1, d_eps_2 = invPartiel(H_elastique(linéaire = False,p = PQ[-1][0])) @ np.array([d_eps_1, d_sig_3])

            Sigma.append(Sigma[-1] + d_sig)
            Epsilon.append(Epsilon[-1] + d_eps)
            PQ.append(Mp@Sigma[-1])
            Epsilon_vp.append(Me@Epsilon[-1])
            F.append(f(PQ[-1], p0))

        else:

            d_sig_1, d_eps_2 = invPartiel(H_plastique(PQ[-1][0]), p0) @ np.array([d_eps_1, d_sig_3])



            Sigma.append(Sigma[-1] + np.array([d_sig_1, d_sig_3]))
            Epsilon.append(Epsilon[-1] + np.array([d_eps_1, d_eps_3]))
            PQ.append(Mp@Sigma[-1])
            Epsilon_vp.append(Me@Epsilon[-1])

            p0 = p0_uptudate(PQ[-1])
            F.append(f(PQ[-1], p0))
    T.append(len(Sigma))
    return np.array(Sigma), np.array(Epsilon), np.array(PQ), np.array(Epsilon_vp), F,p0


def print_graphe_part_2(Sigma, Epsilon, PQ, Epsilon_vp,title,F,p01):

    surfarce_p0 = display_surface_de_charge(p01)

    figure = plt.figure(figsize = (20,10)) # pour afficher les 4 courbes en même temps pour mieux comparer les différentes méthodes
    plt.figure(1)
    plt.suptitle(title)
    plt.subplots_adjust(wspace=0.3, hspace=.3)
    
    plt.subplot(2, 4, 1)
    plt.plot(PQ[1:,0],PQ[1:,0],color='#9F00EC')
    plt.title(' $\sigma_3$ en fonction de $\sigma_1$')
    plt.ylabel('$\sigma_3$')
    plt.xlabel('$\sigma_1$')
    
    plt.subplot(2, 4, 2)
    plt.plot(PQ[1:,0],PQ[1:,0],color='green')
    plt.plot(surfarce_p0[1:,0],surfarce_p0[1:,0],color='red')
    plt.title("q en fonction de p")
    plt.ylabel('q')
    plt.xlabel('p')
    
    plt.subplot(2, 4, 3)
    plt.plot(PQ[1:,0],indice_des_vides(Epsilon_vp[1:,0]),color='blue')
    plt.title("e en fonction de p")
    plt.ylabel('e')
    plt.xlabel('p')
    
    plt.subplot(2, 4, 4)
    plt.plot(np.log(PQ[1:,0]),indice_des_vides(Epsilon_vp[1:,0]),color='red')
    plt.title("e en fonction de ln(p) ")
    plt.ylabel('e')
    plt.xlabel('ln(p)')
    
    plt.subplot(2, 4, 5)
    plt.plot(Epsilon[1:,0],Epsilon_vp[1:,0],color='red')
    plt.title("$\epsilon_v$ en fonction de $\epsilon_1$")
    plt.ylabel('$\epsilon_3$')
    plt.xlabel('$\epsilon_1$')
    
    plt.subplot(2, 4, 6)
    plt.plot(Epsilon[1:,0],PQ[1:,0],color='green')
    plt.title(" q en fonction de $\epsilon_1$")
    plt.ylabel('q')
    plt.xlabel('$\epsilon_1$')
    

    plt.subplot(2, 4, 7)
    plt.plot(Epsilon[1:,],F[1:],color='green')
    plt.title(" f en fonction de $\epsilon_1$")
    plt.ylabel('f')
    plt.xlabel('$\epsilon_1$')

    #plt.savefig(title+'.png')
    plt.show()

    
    return None


def partie2():

    Sigma, Epsilon, PQ, Epsilon_vp, F, p0f = essai_oedométrique()
    print_graphe_part_2(Sigma, Epsilon, PQ, Epsilon_vp,"Essai Oedometrique",F,p0f)
    plt.show()

    Sigma, Epsilon, PQ, Epsilon_vp, F, p0f= essai_isotrope_draine()
    print_graphe_part_2(Sigma, Epsilon, PQ, Epsilon_vp,"Essai isotropique drainé",F,p0f)
    plt.show()

    Sigma, Epsilon, PQ, Epsilon_vp, F,p0f= essai_triaxial_draine()
    print_graphe_part_2(Sigma, Epsilon, PQ, Epsilon_vp,"Essai triaxial drainé",F,p0f)
    plt.show()

    Sigma, Epsilon, PQ, Epsilon_vp, F,p0f= essai_triaxial_draine_surconso()
    print_graphe_part_2(Sigma, Epsilon, PQ, Epsilon_vp,"Essai triaxial drainé très consolidé",F,p0f)
    plt.show()




partie2()