#calculate density

import math as m
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad

#hbar and c =1
Z = 20
M = 0.00939272 #mass nucleon in 10^2GeV
A = 40 #mass number Argon
Omega = (45*(A**(-1./3.)))-(25*(A**(-2./3.)))*10**(-5) #parametrization for hbar*Omega for HO in 10^2GeV
v = M*Omega 
v_mom = 1/(Omega*M) 

#quantum numbers
n = 2
l = 2

#returns unnormalized radial wavefunction calculated with constant a 
def Rad_wavefunc(k, a, n, l):
	t = a*(k**2)
	return (k**l)*sp.special.eval_genlaguerre(n, l+0.5, t)*np.exp(-(a/2.)*k**2)

def norm(a, n, l):
	return np.sqrt((2*m.factorial(n)*a**(l+(3./2.)))/m.gamma(n+l+(3./2.)))

def one_b_mom_dist(k):	
	i = 0
	j = 0
	result = 0
	while (i < Z):
		while (j < n+1):
			result = result + (norm(v_mom, i, j)**2)*Rad_wavefunc(k, v_mom, i, j)**2
			j = j+1 
		i = i+1
	return result

f = open("oneb_momentum_dist.dat", "w")	
i = 0	
while (i < 3000):
	txt = str(one_b_mom_dist(0.001*i)) + '\t' + str(0.001*i)
	f.write( txt  ) 
	f.write('\n')
	i = i+1	