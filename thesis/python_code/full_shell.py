#two-body momentum for full shell nucleus
import math
from moshinsky import Moshinsky, clebsch_gordan
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad


nuc = ["He","O","Ca40"]
A   = [4,16,40]
Z   = [2,8,20]
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
     return (r**2)*Rad_wavefunc(r, v(M_w,0), n, l)*norm(v(M_w,0), n, l)*sph_Bessel(l, k*r)
     
#integral thesis M. Vanhalst p.75 eq.16
def wave_part(n, l, M_w,k): 
     result, err = sp.integrate.quad(integr,0, 100, args=(n,l, M_w, k))
     return result
     
def coefficient(S,M_S,T,M_T, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N, L, M_L):
     result = 0
     for lambd in xrange(abs(l_1 - l_2), l_1 + l_2 +1):
          for M_lambd in xrange(-lambd, lambd +1):
               result = result + clebsch_gordan(l_1, l_2, lambd, m_1, m_2, M_lambd)*Moshinsky(n,l,N,L,n_1,l_1,n_2,l_2,lambd)*clebsch_gordan(l, L, lambd, m_l, M_L, M_lambd)
     return (1./math.sqrt(2))*(1-(-1)**(S+T+l))*clebsch_gordan(0.5, 0.5, S, s_1, s_2, M_S)*clebsch_gordan(0.5, 0.5, T, t_1, t_2, M_T)*result
     
def twobody(P,n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, i):
     result = 0
     for N in xrange(0, 2*n_1 + l_1 +2*n_2 + l_2 + 1):
          for L in xrange(0, 2*n_1 + l_1 +2*n_2 + l_2 + 1):
               a = Rad_wavefunc(P, v_mom(M_average,i), N, L)*norm(v_mom(M_average,i), N, L)
               for N_ in xrange(0, 2*n_1 + l_1 +2*n_2 + l_2 + 1):
                    b = Rad_wavefunc(P, v_mom(M_average,i), N_, L)*norm(v_mom(M_average,i), N_, L)
                    for n in xrange(0, 2*n_1 + l_1 +2*n_2 + l_2 + 1):
                         for l in xrange(0, 2*n_1 + l_1 +2*n_2 + l_2 + 1):
                              if(2*n_1 + l_1 +2*n_2 + l_2  == 2*n + l + 2*N + L and 2*n_1 + l_1 +2*n_2 + l_2 == 2*n + l + 2*N_ + L):
                                   for m_l in xrange(-l, l+1):
                                        for M_L in xrange(-L,L+1):
                                             for S in xrange(0,2):
                                                  for T in xrange(0,2):
                                                       result = result + coefficient(S, s_1 + s_2, T, t_1 + t_2, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N, L, M_L)*coefficient(S, s_1 + s_2, T, t_1 + t_2, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N_, L, M_L)*a*b
     return (1.0/(A[i]*(A[i]-1)))*result

shell = [0,1,2]
     
def opvulling(P,i):
     result = 0
     for N_1 in xrange(0,shell[i]+1):
          if(N_1 % 2 == 0):
               start_1 = 0
          else:
               start_1 = 1
          for l_1 in xrange(start_1,shell[i]+1,2):
               n_1 = int((N_1 - l_1)/2)
               for m_1 in xrange(-l_1, l_1 +1):
                    for N_2 in xrange(0,shell[i]+1):
                         if(N_2 % 2 == 0):
                              start_2 = 0
                         else:
                              start_2 = 1
                         for l_2 in xrange(start_2,shell[i]+1, 2):
                              n_2 = int((N_2 - l_2)/2)
                              print "n_1 = " + str(n_1)
                              print "n_2 = " + str(n_2)
                              for m_2 in xrange(-l_2, l_2 +1):
                                   for s_1 in xrange(0,2):
                                        for s_2 in xrange(0,2):
                                             for t_1 in xrange(0,2):
                                                  for t_2 in xrange(0,2):
                                                       if(n_1 != n_2 or l_1 != l_2 or m_1 != m_2 or s_1 != s_2 or t_1 != t_2):
                                                            result = result + twobody(P,n_1, l_1, m_1, s_1-0.5, t_1-0.5, n_2, l_2, m_2, s_2-0.5, t_2-0.5, i)
     return result
     

def format_(value):
    return "%.8f" % value

distr_array = []
step = 0.1
upperlimit = 26
for j in xrange(0,2):
     f = open("{:s}.txt".format(nuc[j]), "w")
     f.write("# mass number (A) = " + str(A[j]) +  '\n')
     f.write("# number of protons (Z) = " + str(Z[j]) +'\n')
     f.write("# k (1/fm)" + '\t' + "total two_body (fm^3)"  +  '\n')
     i = 0
     while (i < upperlimit):
          print i
          c = 0
          c = opvulling(step*i, j)
          distr_array.append(c*((step*i)**2))
          txt = str(step*i) + '\t' + str(format_(c)) 
          f.write(txt)
          f.write('\n')
          i = i+1
     nor = sp.integrate.trapz(distr_array,x=None,dx = step)
     f.write("#normalisation = " + str(nor) +'\n') 
     f.close()

                                                       
                                   
                              
                         
                         
                         
                

     
     
     
     
     
     