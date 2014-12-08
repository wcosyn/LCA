import math as ma
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad


#quantum numbers
A = int(input("mass number: "))
Z = int(input("number of protons: "))
N = A - Z #number of neutrons
n = int(input("n: "))
l = int(input("l: "))
m = int(input("m: "))
hbarc = 197.327 #(MeV*fm)
M = 939 #mass nucleon (MeV) 
Omega = (45*(A**(-1./3.)))-(25*(A**(-2./3.))) #parametrization for hbar*Omega for HO (MeV)
v = M*Omega / (hbarc**2)
v_mom = (hbarc**2) / (Omega*M)

def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)

#returns unnormalized radial wavefunction calculated with constant a 
def Rad_wavefunc(k, a, n, l):
	t = a*(k**2)
	return (k**l)*sp.special.eval_genlaguerre(n, l+0.5, t)*np.exp(-(a*0.5)*k**2)

def norm(a, n, l):
	return ((2*factorial(n)*a**(l+(3./2.)))/ma.gamma(n+l+(3./2.)))**(0.5)

f = open("data.txt", "w")
i = 0
while (i < 3000):
     txt = str(norm(v_mom, n, l)*Rad_wavefunc(i*10**(-3), v_mom, n, l)) + '\t' + str(i*10**(-3))
     f.write( txt  ) 
     f.write('\n')
     i = i+ 1
f.close()


def ft_integrand(r, n, l, k):
     return (4*np.pi/(2*np.pi)**(3./2.))*norm(v, n, l)*Rad_wavefunc(r, v, n, l)*(r**2)*(np.pi/(2*k*r))**(0.5)*sp.special.jn(l+0.5, k*r)


def FT(n, l, k):
     ans, err = quad(ft_integrand, 0,10, args=(n, l, k))
     return ans

f = open("data2.txt", "w")	
i = 1	
while (i<3000):
          txt = str(FT(n, l, i*10**(-3))) + '\t' + str(i*10**(-3))
          f.write( txt  ) 
          f.write('\n')
          i = i + 1
f.close()



