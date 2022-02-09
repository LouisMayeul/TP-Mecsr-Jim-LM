import numpy as np
import Elastique_LM as elast


p_ref = 1 * 1e3
e0 = 0.8
p_0 = 100 * 1e3
lamdba = 0.2
kappa = 0.01
M = 1.0
G = 2 * 1e6



def M_plastique(p,q):
    a = (2*p - p_0)/(p*p_0*(1 + e0)/(lamdba - kappa)) + kappa/((1 + e0)*p)
    b = 2*q/(M**2*p*p_0*(1 + e0)/(lamdba - kappa))
    c = 2*q*(2*p - p_0)/(M**2*p*p_0*(1 + e0)/(lamdba - kappa)*(2*p - p_0))
    d = 4*q**2/(M**4*p*p_0*(1 + e0)/(lamdba - kappa)*(2*p - p_0)) + 1/(3*G)
    return np.array([[a, b], [c, d]])

def f(p,q):
    return q**2 - M**2*p*(p_0 - p)



"""
########################################################
###################### Question 1 ######################
########################################################
"""


sigma_1 = 10 * 1e3
sigma_2 = sigma_1/2
sigma_3 = sigma_1/2


def Essai_oedometrique(pas_eps = 1e-5, linéaire = True):
    eps_max = 0.2
    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([sigma_1, sigma_2])]
    PQ      = [elast.Mp @ Sigma[-1]]
    Epsilon_vp = [np.array([0,0])]
    
    preconsolidation = p_0

    while Epsilon[-1][0]<eps_max:
        #calcul
        eps = np.array([Epsilon[-1][0] + pas_eps, Epsilon[-1][1]])
        
        if Sigma[-1][0] < preconsolidation :
            sig = np.linalg.inv(elast.H(linéaire = True, p = PQ[-1][0])) @ eps

        else :
            M = M_plastique(p = PQ[-1][0], q = PQ[-1][1])
            sig = np.linalg.inv(M) @ eps
            preconsolidation = sig[0]
        
        
        #Stock
        Sigma.append(sig)
        Epsilon.append(eps)
        PQ.append(elast.Mp@sig)
        Epsilon_vp.append(elast.Me@eps)

    return np.array(Sigma), np.array(Epsilon), np.array(PQ), np.array(Epsilon_vp)

Sigma, Epsilon, PQ, Epsilon_vp = Essai_oedometrique()
elast.print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,"Essai oedometrique - Deformations")


"""
########################################################
###################### Question 2 ######################
########################################################
"""


sigma_1 = 10 * 1e3
sigma_2 = sigma_1
sigma_3 = sigma_1


def Essai_isotrope_draine(pas_sig = 1e3, linéaire = True):
    sig_max = 250 * 1e3
    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([sigma_1, sigma_2])]
    PQ      = [elast.Mp @ Sigma[-1]]
    Epsilon_vp = [np.array([0,0])]
    
    preconsolidation = p_0

    while Sigma[-1][0]<sig_max:
        #calcul
        sig = np.array([Sigma[-1][0] + pas_sig, Sigma[-1][1] + pas_sig])
        
        if Sigma[-1][0] < preconsolidation :
            eps = elast.H(linéaire = linéaire,p = PQ[-1][0]) @ sig 

        else :
            M = M_plastique(p = PQ[-1][0], q = PQ[-1][1])
            eps = M @ eps
            preconsolidation = sig[0]
    
        
        #Stock
        Sigma.append(sig)
        Epsilon.append(eps)
        PQ.append(elast.Mp@sig)
        Epsilon_vp.append(elast.Me@eps)

    return np.array(Sigma), np.array(Epsilon), np.array(PQ), np.array(Epsilon_vp)

Sigma, Epsilon, PQ, Epsilon_vp = Essai_isotrope_draine()
elast.print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,"Essai isotrope draine - Contraintes")



"""
########################################################
###################### Question 3 ######################
########################################################
"""


sigma_1 = 10 * 1e3
sigma_2 = sigma_1
sigma_3 = sigma_1


def Essai_triaxial_draine(pas_sig = 1e3, linéaire = True):
    Epsilon = [np.array([0,0])]
    Sigma   = [np.array([sigma_1, sigma_2])]
    PQ      = [elast.Mp @ Sigma[-1]]
    Epsilon_vp = [np.array([0,0])]
    
    preconsolidation = p_0

    """ Preconsalidation """
    while Sigma[-1][0] < 250 * 1e3:
        #calcul
        sig = np.array([Sigma[-1][0] + pas_sig, Sigma[-1][1] + pas_sig])
        
        if Sigma[-1][0] <= preconsolidation :
            eps = elast.H(linéaire = linéaire,p = PQ[-1][0]) @ sig 
        else :
            M = M_plastique(p = PQ[-1][0], q = PQ[-1][1])
            eps = M @ eps
            preconsolidation = sig[0]
        
        print("sig = ", sig[0]/1000, " et preconsolidation = ", preconsolidation/1000)
        
        #Stock
        Sigma.append(sig)
        Epsilon.append(eps)
        PQ.append(elast.Mp@sig)
        Epsilon_vp.append(elast.Me@eps)
    
    """" décharge isotrope """
    while Sigma[-1][0] > 200 * 1e3:
        #calcul
        sig = np.array([Sigma[-1][0] - pas_sig, Sigma[-1][1] - pas_sig])
        
        if Sigma[-1][0] <= preconsolidation :
            eps = elast.H(linéaire = linéaire,p = PQ[-1][0]) @ sig 
        else :
            M = M_plastique(p = PQ[-1][0], q = PQ[-1][1])
            eps = M @ eps
            preconsolidation = sig[0]
        
        print("sig = ", sig[0]/1000, " et preconsolidation = ", preconsolidation/1000)
        
        #Stock
        Sigma.append(sig)
        Epsilon.append(eps)
        PQ.append(elast.Mp@sig)
        Epsilon_vp.append(elast.Me@eps)
    
    """ Cisaillement """
    while f(PQ[-1][0], PQ[-1][1]) <= 0: # remplacer par la condition de rupture
        #calcul
        sig = np.array([Sigma[-1][0] + pas_sig, Sigma[-1][1]])
        
        if Sigma[-1][0] <= preconsolidation :
            eps = elast.H(linéaire = linéaire,p = PQ[-1][0]) @ sig 
        else :
            M = M_plastique(p = PQ[-1][0], q = PQ[-1][1])
            eps = M @ eps
            preconsolidation = sig[0]
        
        print("sig = ", sig[0]/1000, " et preconsolidation = ", preconsolidation/1000)
        
        #Stock
        Sigma.append(sig)
        Epsilon.append(eps)
        PQ.append(elast.Mp@sig)
        Epsilon_vp.append(elast.Me@eps)

    return np.array(Sigma), np.array(Epsilon), np.array(PQ), np.array(Epsilon_vp)

Sigma, Epsilon, PQ, Epsilon_vp = Essai_triaxial_draine()
elast.print_graphe(Sigma, Epsilon, PQ, Epsilon_vp,"Essai triaxial draine - Contraintes")




"""
########################################################
###################### Question 4 ######################
########################################################
"""

