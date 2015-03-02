#two-body momentum for Helium (state 1s_1/2) (Clebsch_gordan = 1, Moshinsky = 1)

import math
from moshinsky import Moshinsky, clebsch_gordan
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad


nuc = ["He","Be","C","O","Al","Ca40","Ca48","Fe","Ag","Ar","Pb"]
A   = [4,9,12,16,27,40,48,56,109,40,208]
Z   = [2,4,6,8,13,20,20,26,47,18,82]
N = [m - n for m,n in zip(A,Z)] #number of neutrons

Omega = []
for i in range(len(A)):
     Omega.append((45*(A[i]**(-1./3.)))-(25*(A[i]**(-2./3.))))
     
hbarc = 197.327 #(MeV*fm)
M_n = 939.565 #mass neutron (MeV) [Particle Data Group]
M_p = 938.272 #mass proton (MeV) [Particle Data Group]
M_average = (M_n + M_p)/2.



#returns HO constant \nu for different mass
def v(M, index):
     return M*Omega[index] / (hbarc**2)
     
#returns HO constant \nu in momentum space for different mass
def v_mom(M, index):
     return (hbarc**2) / (Omega[index]*M) 

#returns unnormalized radial wavefunction for HO calculated with constant a 
def Rad_wavefunc(k, a, n, l):
	t = a*(k**2)
	return (k**l)*sp.special.eval_genlaguerre(n, l+0.5, t)*np.exp(-(a/2.)*(k**2))
     
#norm radial wave function
def norm(a, n, l):
	return np.sqrt((2*math.factorial(n)*(a**(l+(3./2.))))/math.gamma(n+l+(3./2.)))

#returns spherical bessel function of order n in point z     
def sph_Bessel(n,z):
     x = sp.special.sph_jn(n,z)
     return x[0][n]

#integrand for integrals thesis M. Vanhalst p.75 eq.16
def integr(r,n,l, M_w,k):
     return (r**2)*Rad_wavefunc(r*np.sqrt(2.), v(M_w,0), n, l)*norm(v(M_w,0), n, l)*sph_Bessel(l, k*r)
     
#integral thesis M. Vanhalst p.75 eq.16
def wave_part(n, l, M_w,k): 
     result, err = sp.integrate.quad(integr,0, 100, args=(n,l, M_w, k))
     return result

def coefficient(S,M_S,T,M_T, s_1, s_2, t_1, t_2):
     return (1./math.sqrt(2))*(1-(-1)**(S+T))*clebsch_gordan(0.5, 0.5, S, s_1, s_2, M_S)*clebsch_gordan(0.5, 0.5, T, t_1, t_2, M_T)
     
def two_body(P):
     result = 0
     for ms_1 in xrange(0,2):
          for ms_2 in xrange(0,2):
               for mt_1 in xrange(0,2):
                    for mt_2 in xrange(0,2):
                         for S in xrange(0,2):
                              for M_S in xrange(-S,S+1):
                                   for T in xrange(0,2):
                                        for M_T in xrange(-T,T+1):
                                             if(ms_1 != ms_2 or mt_1 != mt_2):
                                                  result = result + (coefficient(S,M_S,T,M_T,ms_1-0.5,ms_2-0.5,mt_1-0.5,mt_2-0.5)**2)*(wave_part(0,0, M_average,P)**2)
     return (2.0/math.pi)*(1.0/12.0)*result

def format_(value):
    return "%.8f" % value

distr_array = []
step = 0.1
upperlimit = 50
f = open("helium.txt", "w")
f.write("# mass number (A) = " + str(A[0]) +  '\n')
f.write("# number of protons (Z) = " + str(Z[0]) +'\n')
f.write("# k (1/fm)" + '\t' + "total two_body (fm^3)"  +  '\n')
i = 0
while (i < upperlimit):
     c = 0
     c = two_body(step*i)
     distr_array.append(c*(step*i)**2)
     txt = str(step*i) + '\t' + str(format_(c)) 
     f.write(txt)
     f.write('\n')
     i = i+1
nor = sp.integrate.trapz(distr_array,x=None,dx = step)
f.write("#normalisation = " + str(nor) +'\n') 
f.close()

     
     
                                   
     