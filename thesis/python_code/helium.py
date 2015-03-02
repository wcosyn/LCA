#calculate two-body total momentum (\vec{P} = \vec{k_1} + \vec{k_2}) distribution 

import math
from moshinsky import clebsch_gordan, Moshinsky 
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
	return (k**l)*sp.special.eval_genlaguerre(n, l+0.5, t)*np.exp(-(a/2.)*k**2)
     
#norm radial wave function
def norm(a, n, l):
	return np.sqrt((2*math.factorial(n)*a**(l+(3./2.)))/math.gamma(n+l+(3./2.)))
     
def integrand_norm(k,n,l):
     return (Rad_wavefunc(k, v_mom(M_average,0), n, l)*norm(v_mom(M_average,0), n, l)*k)**2

     
def normalis(n,l):
     nor, err = sp.integrate.quad(integrand_norm,0,20, args= (n,l))
     return nor
     

print normalis(0,0)    
#returns spherical bessel function of order n in point z     
def sph_Bessel(n,z):
     x = sp.special.sph_jn(n,z)
     return x[0][n]
     
f = open("bessel.txt", "w")
i = 0
while (i < 2001):
     d = sph_Bessel(0,0.01*i)
     txt = str(0.01*i) + '\t' + str((d))  
     f.write(txt)
     f.write('\n')
     i = i+1
f.close()
  
                         for S in xrange(0,2):
                              for M_S in xrange(-S,S+1):
#rechtstreekse input van k-space golffuncties Rad_wavefunc(P, v_mom, N, L)
def twobody_1(k, M_w, i, mt_1, mt_2, n1 , n2):
     result = 0
     for ms_1 in xrange(0, 2):
          for ms_2 in xrange(0, 2):
               for S in xrange(0, 2):
                    for T in xrange(0, 2):
                         if(ms_1 != ms_2 or mt_1 != mt_2):
                              result = result + ((1-(-1)**(S+T))**2)*(clebsch_gordan(0.5,0.5,S,ms_1-0.5,ms_2-0.5,ms_1 + ms_2 - 1)**2)*(clebsch_gordan(0.5,0.5,T,mt_1,mt_2,mt_1 + mt_2)**2)*Rad_wavefunc(k/np.sqrt(2), v_mom(M_w,i), 0, 0)*Rad_wavefunc(k/np.sqrt(2), v_mom(M_w,i), 0, 0)*norm(v_mom(M_w,i), 0, 0)*norm(v_mom(M_w,i), 0, 0)
     return (1/16.)*result
     
#berekening eq 16 p75 thesis M. Vanhalst
def twobody_2(k, M_w, i, mt_1, mt_2, n1 , n2):
     result = 0
     for ms_1 in xrange(0, 2):
          for ms_2 in xrange(0, 2):
               for S in xrange(0, 2):
                    for M_S in xrange(-S, S+1):
                         for T in xrange(0, 2):
                              for M_T in xrange(-T, T+1):
                                   if(ms_1 != ms_2 or mt_1 != mt_2):
                                        result = result + ((1-(-1)**(S+T))**2)*(clebsch_gordan(0.5,0.5,S,ms_1-0.5,ms_2-0.5,M_S)**2)*(clebsch_gordan(0.5,0.5,T,mt_1,mt_2,M_T)**2)*wave_part(0,0, M_w, i,k)*wave_part(0, 0, M_w, i,k)
     return (1./math.pi)*result
     
def opvulling_1(k,Q_1, Q_2, M_w, i, mt_1, mt_2):
     return twobody_1(k, M_w, i, mt_1, mt_2,0, 0)
     
#twobody summation for pp, nn and pn
def distribution_1(k, i):
     return (1.0/(A[i]*(A[i]-1)))*(opvulling_1(k, Z[i], Z[i],M_average, i, 0.5, 0.5) + opvulling_1(k,N[i], N[i], M_average, i, -0.5, -0.5) + opvulling_1(k, Z[i], N[i],M_average, i, 0.5, -0.5) + opvulling_1(k, N[i], Z[i],M_average, i, -0.5, 0.5))
     

def opvulling_2(k,Q_1, Q_2, M_w, i, mt_1, mt_2):
     return twobody_2(k, M_w, i, mt_1, mt_2,0, 0)
     
#twobody summation for pp, nn and pn
def distribution_2(k, i):
     return (1.0/(A[i]*(A[i]-1)))*(opvulling_2(k, Z[i], Z[i],M_average, i, -0.5, -0.5) + opvulling_2(k,N[i], N[i], M_average, i, 0.5, 0.5) + opvulling_2(k, Z[i], N[i],M_average, i, -0.5, 0.5) + opvulling_2(k, N[i], Z[i],M_average, i, 0.5, -0.5))
     
def integrand(k,i):
     return distribution_1(k,i)*(k**2)
#helpfunction: round values
def format_(value):
    return "%.8f" % value

nuclides = [0]
distr_array_1 = []
distr_array_2 = []
step = 0.1
upperlimit = 25
for j in nuclides:
     f = open("{:s}_mf_twobody.txt".format(nuc[j]), "w")
     f.write("# mass number (A) = " + str(A[j]) +  '\n')
     f.write("# number of protons (Z) = " + str(Z[j]) +'\n')
     f.write("# k (1/fm)" + '\t' + "mf_zelf (fm^3)" + '\t' + "mf_maarten (fm^3)" + '\n')
     print str(A[j])
     i = 0
     while (i < upperlimit):
          c = distribution_1(step*i, j)
          d = distribution_2(step*i, j)
          distr_array_1.append(c*(step*i)**2)
          distr_array_2.append(d*(step*i)**2)
          txt = str(step*i) + '\t' + str(format_(c)) + '\t' + str(format_(d))  + '\t' + str(format_(c/d))
          f.write(txt)
          f.write('\n')
          i = i+1
     nor_1 = sp.integrate.trapz(distr_array_1,x=None,dx = step)
     nor_2 = sp.integrate.trapz(distr_array_2,x=None,dx = step)
     nor = sp.integrate.quad(integrand, 0, 5, args=(j) )
     f.write("#normalisation = " + str(nor_1) + '\t' + str(nor_2) + '\t' + str(nor) +'\n') 
     f.close()

   
     
     