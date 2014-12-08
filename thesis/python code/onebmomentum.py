#calculate one-body momentum distribution

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
     

	
def one_b_mom_dist(k, Q, M):
     f = 0	
     n = 0
     counter = 0
     result = 0
     l = 0
     d = 0
     N = 0 
     while (counter < Q):
          if N % 2 == 0:
               l = 0
          else:
               l = 1
          while (l <= N and counter < Q):
               d = 2*((2*l)+1)
               n = (N-l)/2
               if counter + d < Q+1:
                    f = 1
               else:
                    f = (float(Q - counter))/float(d)
               result = result + float(d)*float(f)*(norm(v_mom(M), n, l)*Rad_wavefunc(k, v_mom(M), n, l))**2
               counter = int(counter + d*f)
               l = l + 2
     return result
                    
                             
def distribution(k):
     return (one_b_mom_dist(k, N, M_n)+ one_b_mom_dist(k, Z, M_p))/A


f = open("oneb_momentum_dist.txt", "w")	
i = 0	
while (i < 40000):
     txt = str(distribution(0.0001*i)) + '\t' + str(0.0001*i)
     f.write( txt  )
     f.write('\n')
     i = i+1
f.close()
