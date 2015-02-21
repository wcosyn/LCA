#calculate two-body momentum distribution

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

def norm(a, n, l):
	return np.sqrt((2*math.factorial(n)*a**(l+(3./2.)))/math.gamma(n+l+(3./2.)))
     
    
#opvulling van Wood-Saxon
#degeneracy
degen   = np.array([2,4,2,6,2,4,8,4,6,2,10,8,6,4,2,12,10,8,6,4,2,14,10,6,12,8,2,4])

#orbital Quantum number l
orbital = np.array([0,1,1,2,0,2,3,1,3,1,4,4,2,2,0,5,5,3,3,1,1,6,4,2,6,4,0,2])

#radial Quantum number n
radial = np.array([0,0,0,0,1,0,0,1,0,1,0,0,1,1,2,0,0,1,1,2,2,0,1,2,0,1,3,2])



#two-body distribution (k is relative momentum) for Q_1 particles in one HO well and Q_2 particles in another HO well ( for example neutrons and protons are in different wells) (M_w: mass particles) (mt_1,mt_2: isospin projections)
def twobody(k, M_w, i, mt_1, mt_2, n1 , n2):
     result = 0
     for m1 in xrange(-orbital[n1],orbital[n1] + 1):
          for m2 in xrange(-orbital[n2],orbital[n2] + 1):
               for ms_1 in xrange(0, 2):
                    for ms_2 in xrange(0, 2):
                         for n in xrange(0, 2*radial[n1] + orbital[n1] + 2*radial[n2] + orbital[n2] + 1):
                              for n_ in xrange(0, 2*radial[n1] + orbital[n1] + 2*radial[n2] + orbital[n2] + 1):
                                   for l in xrange(0, orbital[n1] + orbital[n2] + 1):
                                        for m_l in xrange(-l,l+1):
                                             for N in xrange(0, 2*radial[n1] + orbital[n1] + 2*radial[n2] + orbital[n2] + 1):
                                                  for L in xrange(0, orbital[n1] + orbital[n2] + 1):
                                                       for M_L in xrange(-L,L+1):
                                                            for S in xrange(0, 2):
                                                                 for M_S in xrange(-S, S+1):
                                                                      for T in xrange(0, 2):
                                                                           for M_T in xrange(-T, T+1):
                                                                                for lam_p in xrange(int(math.fabs(l-L)), l+L+1):
                                                                                     for M_lam_p in xrange(-lam_p,lam_p+1):
                                                                                          for lam in xrange(int(math.fabs(l-L)), l+L+1):
                                                                                               for M_lam in xrange(-lam,lam+1):
                                                                                                    if(2*radial[n1] + orbital[n1] + 2*radial[n2] + orbital[n2] == 2*n+l+2*N+L and 2*radial[n1] + orbital[n1] + 2*radial[n2] + orbital[n2] == 2*n_ +l+2*N+L ):
                                                                                                         result = result + ((1-(-1)*(l+S+T))**2)*(clebsch_gordan(0.5,0.5,S,ms_1-0.5,ms_2-0.5,M_S)**2)*(clebsch_gordan(0.5,0.5,T,mt_1,mt_2,M_T)**2)*(clebsch_gordan(orbital[n1],orbital[n2],lam,m1,m2,M_lam)*clebsch_gordan(orbital[n1],orbital[n2],lam_p,m1,m2,M_lam_p))*(clebsch_gordan(l,L,lam,m_l,M_L,M_lam)*clebsch_gordan(l,L,lam_p,m_l,M_L,M_lam_p))*Moshinsky(n,l,N,L,radial[n1],orbital[n1],radial[n2], orbital[n2],lam)*Moshinsky(n_,l,N,L,radial[n1],orbital[n1],radial[n2], orbital[n2],lam_p)*Rad_wavefunc(k, v_mom(M_w,i), n, l)*Rad_wavefunc(k, v_mom(M_w,i), n_, l)*norm(v_mom(M_w,i), n, l)*norm(v_mom(M_w,i), n_, l)
     return (0.5)*result
     
 
def opvulling(k,Q_1, Q_2, M_w, i, mt_1, mt_2):
     z = 0
     j = 0
     f = 0
     b = 0
     result = 0
     while(b < Q_2):
          b = b + degen[i]
          z = z + 1
     while(f < Q_2):
          f = f + degen[j]
          j = j + 1
     for n1 in xrange(0,z):
          for n2 in xrange(0,j):
               result = result + twobody(k, M_w, i, mt_1, mt_2,n1, n2)
     return result
#twobody summation for pp, nn and pn
def distribution(k, i):
     return (1.0/(A[i]*(A[i]-1)))*(opvulling(k, Z[i], Z[i],M_average, i, 0.5, 0.5) + opvulling(k,N[i], N[i], M_average, i, -0.5, -0.5) + opvulling(k, Z[i], N[i],M_average, i, 0.5, -0.5))
     
  
#helpfunction: round values     
def format(value):
    return "%.8f" % value

nuclides = [0,2]
distr_array = []
for j in nuclides:
     f = open("{:s}_mf_twobody.txt".format(nuc[j]), "w")
     f.write("# mass number (A) = " + str(A[j]) +  '\n')
     f.write("# number of protons (Z) = " + str(Z[j]) +'\n')
     f.write("# k (1/fm)" + '\t' + "mf (fm^3)" + '\n')
     print str(A[j])
     i = 0
     while (i < 25):
          print i
          c = distribution(0.1*i, j)
          distr_array.append(c*(0.1*i)**2)
          txt = str(0.1*i) + '\t' + str(format(c))
          f.write(txt)
          f.write('\n')
          i = i+1
     nor = sp.integrate.trapz(distr_array,x=None,dx = 0.1)
     f.write("#normalisation = " + str(nor) + '\n') 
     f.close()

         
     
     
     
     