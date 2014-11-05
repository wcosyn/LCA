#calculate density

import math as m
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad

hbarc = 197.327 #(MeV*fm)
M = 939.272 #mass nucleon (MeV)
A = 40 #mass number Argon
Z = 3
Omega = (45*(A**(-1./3.)))-(25*(A**(-2./3.))) #parametrization for hbar*Omega for HO (MeV)
v = M*Omega / hbarc
v_mom = hbarc / (Omega*M) 

#returns unnormalized radial wavefunction calculated with constant a 
def Rad_wavefunc(k, a, n, l):
	t = a*(k**2)
	return (k**l)*sp.special.eval_genlaguerre(n, l+0.5, t)*np.exp(-(a/2.)*k**2)

def norm(a, n, l):
	return np.sqrt((2*m.factorial(n)*a**(l+(3./2.)))/m.gamma(n+l+(3./2.)))


vf = open("mom_wave_function.txt", "w")	
i = 0	
while (i < 60000):
	txt = str(norm(v_mom, 2, 1)*Rad_wavefunc(0.001*i, v_mom, 2,1)) + '\t' + str(0.001*i)
	vf.write( txt  ) 
	vf.write('\n')
	i = i+1
	
def one_b_mom_dist(k):	
	i = 0
	j = 0
	result = 0
	while (i < Z+1):
		while (j < i+1):
			result = result + (norm(v_mom, i, j)**2)*Rad_wavefunc(k, v_mom, i, j)**2
			j = j+1 
		i = i+1
	return result

f = open("oneb_momentum_dist.txt", "w")	
i = 0	
while (i < 60000):
	txt = str(one_b_mom_dist(0.001*i)) + '\t' + str(0.001*i)
	f.write( txt  ) 
	f.write('\n')
	i = i+1	
