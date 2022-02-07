import Elastique_LM

M = 1
p0 = 1 
lamdba = 1
k=1
e0 = 1


def M_plastique(p,q):
    a = M**2*(2*p-p0)/(M**2*p*p0*(1+e0)/(lamdba -k)) + k(1+e0)*p
    b = 2*q/(M**2*p*p0*(1+e0)/(lamdba-k))
    c = 2*q*M**2*(2*p-p0)/(M**4*p*p0*(1+e0)/(lamdba-k)*(2*p-p0))
    d = 1