#momentum wave functions
import math as ma
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad


#quantum numbers
A = 40
Z = 18
N = A - Z #number of neutrons
n = int(input("n: "))
l = int(input("l: "))
space = bool(input("Space or Momentum wave function (True or False):" ))
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

f = open("wave_function.txt", "w")
i = 0
t = 0
if (space == True):
     t = v
elif (space == False):
     t = v_mom
     
while (i < 8000):
     txt = str(norm(t, n, l)*Rad_wavefunc(i*10**(-3), t, n, l)) + '\t' + str(i*10**(-3))
     f.write( txt  ) 
     f.write('\n')
     i = i+ 1
f.close()