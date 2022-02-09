import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi, sqrt

# Loi de Hooke

v = 0.3
E = 9000
n = 10
G = 2 * 10**6 #*Pa
e0 = 0.80
k = 0.01



H = 1 / E * np.array([[1, -2 * v], [-v, 1 - v]])
Me = np.array([[1, 2], [2 / 3, -2 / 3]])
Mp = np.array([[1 / 3, 2 / 3], [1, -1]])

def inv_part(H):
    a=H[0][0]
    b=H[0][1]
    c=H[1][0]
    d=H[1][1]
    return np.array([[1/a,-b/a],[c/a,d-b*c/a]]) #[[1/a,-b/a],[c/a,d-cb/a]] #inversion partielle

def hooke2(sig10, sig30): #lnp e faux et unités etranges 
    sig1 = sig10
    sig3 = sig30

    H = 1 / E * np.array([[1, -2 * v], [-v, 1 - v]])
    H_inv = inv_part(H)

    eps1 = 0
    eps3 = 0

    epsv=0
    epsq=0

    PQ=Mp@np.array([sig1,sig3])
    p=PQ[0]
    q=PQ[1]

    liste_e = [e0]

    liste_p = [p]
    liste_q = [q]
    liste_sig1 = [sig1]
    liste_sig3 = [sig3]
    liste_eps1 = [eps1]
    liste_epsv = [epsv]


    while sig1 < 500:
        deps1 = n  # déformation axiale et contrainte radiale controlées
        dsig3 = 0

        inconnues=H_inv@np.array([deps1,dsig3])
        dsig1=inconnues[0]
        deps3=inconnues[1]


        dEPSvq = Me@np.array([deps1, deps3])
        depsv = dEPSvq[0]
        depsq = dEPSvq[1]

        dPQ = Mp@np.array([dsig1, dsig3])
        dp = dPQ[0]
        dq = dPQ[1]

        sig1 += dsig1
        sig3 += dsig3
        eps1 += deps1
        eps3 += deps3
        epsq += depsq
        epsv += depsv
        p += dp
        q += dq

        liste_p.append(p)
        liste_q.append(q)
        liste_sig1.append(sig1)
        liste_sig3.append(sig3)
        liste_eps1.append(eps1)
        liste_epsv.append(epsv)

        liste_e.append(-epsv * (1 + e0) + e0)

    liste_p = np.array(liste_p)
    liste_lnp = np.log(liste_p)

    plt.plot(liste_p, liste_q)
    plt.xlabel("p")
    plt.ylabel("q")
    plt.show()
    plt.plot(liste_sig1, liste_sig3)
    plt.xlabel("$\sigma_1$")
    plt.ylabel("$\sigma_3$")
    plt.show()
    plt.plot(liste_eps1, liste_epsv)
    plt.xlabel("$\epsilon_1$")
    plt.ylabel("$\epsilon_v$")
    plt.show()
    plt.plot(liste_eps1, liste_q)
    plt.xlabel("$\epsilon_1$")
    plt.ylabel("q")
    plt.show()
    plt.plot(liste_p, liste_e)
    plt.xlabel("p")
    plt.ylabel("e")
    plt.show()
    plt.plot(liste_lnp, liste_e)
    plt.xlabel("ln(p)")
    plt.ylabel("e")
    plt.show()

    return sig1, sig3, eps1, eps3, epsv, epsq

def nonlineaireb(sig10, sig30):
    sig1 = sig10
    sig3 = sig30

    PQ = np.dot(Mp, np.array([sig1, sig3]))
    p = PQ[0]
    q = PQ[1]

    eps1 = 0
    eps3 = 0
    EPSvq=np.dot(Me,np.array([eps1,eps3]))
    epsv = EPSvq[0]
    epsq = EPSvq[1]

    liste_p = [p]
    liste_q = [q]
    liste_sig1 = [sig1]
    liste_sig3 = [sig3]
    liste_eps1 = [eps1]
    liste_epsv = [epsv]

    liste_e=[e0]


    while sig1 < 500:
        dsig1 = n  # on incrémente le chargement
        dsig3 = 0

        H = np.array([[k / (1 + e0) / p, 0], [0, 1 / 3 / G]])

        M=np.linalg.inv(Me)@H@Mp

        dEPS=M@np.array([dsig1,dsig3])
        deps1 = dEPS[0]
        deps3 = dEPS[1]

        dPQ = Mp@np.array([dsig1, dsig3])
        dp = dPQ[0]
        dq = dPQ[1]

        dEPSvq=Me@dEPS
        depsv = dEPSvq[0]
        depsq = dEPSvq[1]


        sig1 += dsig1
        sig3 += dsig3
        eps1 += deps1
        eps3 += deps3
        epsq += depsq
        epsv += depsv
        p += dp
        q += dq

        liste_p.append(p)

        liste_q.append(q)
        liste_sig1.append(sig1)
        liste_sig3.append(sig3)
        liste_eps1.append(eps1)
        liste_epsv.append(epsv)

        liste_e.append(-epsv*(1+e0)+e0)

    liste_p = np.array(liste_p)
    liste_lnp = np.log(liste_p)

    plt.plot(liste_p, liste_q)
    plt.xlabel("p")
    plt.ylabel("q")
    plt.show()
    plt.plot(liste_sig1, liste_sig3)
    plt.xlabel("$\sigma_1$")
    plt.ylabel("$\sigma_3$")
    plt.show()
    plt.plot(liste_eps1, liste_epsv)
    plt.xlabel("$\epsilon_1$")
    plt.ylabel("$\epsilon_v$")
    plt.show()
    plt.plot(liste_eps1, liste_q)
    plt.xlabel("$\epsilon_1$")
    plt.ylabel("q")
    plt.show()
    plt.plot(liste_p, liste_e)
    plt.xlabel("p")
    plt.ylabel("e")
    plt.show()
    plt.plot(liste_lnp, liste_e)
    plt.xlabel("ln(p)")
    plt.ylabel("e")
    plt.show()

    return sig1, sig3, eps1, eps3, epsv, epsq, p, q