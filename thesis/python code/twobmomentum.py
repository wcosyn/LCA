#calculate two-body momentum distribution

import Wigner as wgr
import math as m
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad


A = int(input("mass number: "))
Z = int(input("number of protons: "))
N = A - Z #number of neutrons

hbarc = 197.327 #(MeV*fm)
M_n = 939.565 #mass neutron (MeV) 
M_p = 938.272 #mass proton (MeV)
Omega = (45*(A**(-1./3.)))-(25*(A**(-2./3.))) #parametrization for hbar*Omega for HO (MeV)

#returns HO constant \nu for different mass
def v(M):
     return M*Omega / (hbarc**2)

def v_mom(M):
     return (hbarc**2) / (Omega*M) 

#returns unnormalized radial wavefunction calculated with constant a 
def Rad_wavefunc(k, a, n, l):
	t = a*(k**2)
	return (k**l)*sp.special.eval_genlaguerre(n, l+0.5, t)*np.exp(-(a/2.)*k**2)

def norm(a, n, l):
	return np.sqrt((2*m.factorial(n)*a**(l+(3./2.)))/m.gamma(n+l+(3./2.)))


def two_body_distr(k,P):
     f1 = 0	
     n2 = 0
     counter1 = 0
     result = 0
     l1 = 0
     d1 = 0
     N1 = 0
     while (counter1 < Q):
          if N1 % 2 == 0:
               l1 = 0
          else:
               l1 = 1
          while (l1 <= N1 and counter1 < Q):
               d1 = 2*((2*l1)+1)
               n1 = (N1-l1)/2
               if counter1 + d1 < Q+1:
                    f1 = 1
               else:
                    f1 = (float(Q - counter1))/float(d1)
               f2 = 0	
               n2 = 0
               counter2 = 0
               l2 = 0
               d2 = 0
               N2 = 0
               while (counter2 < Q):
                    if N2 % 2 == 0:
                         l2 = 0
                    else:
                         l2 = 1
                    while (l2 <= N2 and counter < Q):
                         d2 = 2*((2*l2)+1)
                         n2 = (N2-l2)/2
                         if counter2 + d2 < Q+1:
                              f2 = 1
                         else:
                              f2 = (float(Q - counter2))/float(d2)
                         lamb = m.abs(l1 - l2)
                         while(labm <= l1 + l2):
                              M_lamb = -lamb
                              while(M_lamb <= lamb):
                                   while(n):
                                        while(g):
                                             
                                        n = n + 1
                                   M_lamb = M_lamb + 1
                              lamb = lamb + 1
                         counter2 = int(counter2 + d2*f2)
                         l2 = l2 + 2
               counter1 = int(counter1 + d1*f1)
               l1 = l1 + 2
     return result




def loops():
     while(n <= 10):
          while(N <= 10):
               while(l <= n):
                    while(L <= N)
               n2 = n2 + 1
     n1 = n1 + 1



