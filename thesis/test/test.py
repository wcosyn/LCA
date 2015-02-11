#test for Helium

import math as m
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad


A = int(input("mass number: "))
Z = int(input("number of protons: "))
N = A - Z #number of neutrons

hbarc = 197.327 #(MeV*fm)
M_n = 939.0 #mass neutron (MeV) 
M_p = 938.0 #mass proton (MeV)
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
     

f = open("test.txt", "w")	
i = 0
	
while (i < 50):
     txt = str(2*(Rad_wavefunc(0.1*i,v_mom(M_p), 0, 0)*norm(v_mom(M_p), 0, 0))**2+2*(Rad_wavefunc(0.1*i,v_mom(M_n), 0, 0)*norm(v_mom(M_n), 0, 0))**2) + '\t' + str(0.1*i)
     f.write( txt  )
     f.write('\n')
     i = i+1
f.close()