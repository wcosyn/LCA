import math as m
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad

hbarc = 197.327 #(MeV*fm)
M = 939.272 #mass nucleon (MeV)
A = 40 #mass number Argon
Omega = (45*(A**(-1./3.)))-(25*(A**(-2./3.)))#parametrization for hbar*Omega for HO (MeV)
v = M*Omega / hbarc
v_mom = hbarc / (Omega*M) 

#quantum numbers
n = 0
l = 0

#returns unnormalized radial wavefunction calculated with constant a 
def Rad_wavefunc(r, a, n, l):
	t = a*(r**2)
	return (r**l)*sp.special.eval_genlaguerre(n, l+0.5, t)*np.exp(-(a/2.)*r**2)

f = open("data.txt", "w")	
i = 0
while (i < 5*10**(9)):
	txt = str(Rad_wavefunc(i*10**(-8), v_mom, n, l)) + '\t' + str(i*10**(-8))
	f.write( txt  ) 
	f.write('\n')
	i = i+ (10**(7))

def norm(a, n, l):
	return np.sqrt((2*m.factorial(n)*a**(l+(3./2.)))/m.gamma(n+l+(3./2.)))


def FT_integrand(r, a, n, l, k):
	return Rad_wavefunc(r, a, n, l)*(r**2)*special.sph_jn(l, k*r)[0][-1]


def FT(a, n, l, k):
	return sp.integrate.quad(FT_integrand, 0, 30,limit=500, args=(a, n, l, k))[0]


def compare(a, b, n, l, k):
	if (Rad_wavefunc(k, b, n, l) != 0):
		return (4*np.pi*(1j)**l/(2*np.pi)**(3./2.))*(norm(a, n, l)/norm(b, n, l))*FT(a, n, l, k)/Rad_wavefunc(k, b, n, l)


f = open("data2.dat", "w")	
i = 0	
while (i < 3000):
	txt = str(compare(v, v_mom, n, l, 0.001*i)) + '\t' + str(0.001*i)
	f.write( txt  ) 
	f.write('\n')
	i = i+1
	
			



