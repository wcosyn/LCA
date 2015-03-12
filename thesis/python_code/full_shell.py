#two-body momentum for full shell nucleus
import math
from matrices import clebsch_gordan, clebsch_gordan_spin
import scipy as sp
import numpy as np
from sympy.physics.quantum.cg import CG
from sympy import S
from scipy.integrate import quad
from read import choose_file


nuc = ["He","O"]
A   = [4,16]
Z   = [2,8]
N = [m - n for m,n in zip(A,Z)] #number of neutrons


Omega = []
for i in range(len(A)):
     a = (45*(A[i]**(-1./3.)))-(25*(A[i]**(-2./3.)))
     Omega.append(a)

     
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


#summation over Lambda, M_lambda of clebsch_gordan(l_1, m_1, l_2, m_2, Lambda, M_lambda) * Moshinsky_bracket(n_1, l_1, n_2, l_2, n, l , N , L, Lambda) * clebsch_gordan(l, m_l, L, M_L, Lambda, M_lambda) (see thesis Maarten p74 eq 14)
def CG_MB_CG(n_1, l_1, m_1, n_2, l_2, m_2, n,l , m_l, N, L, M_L):
     result = 0
     z = choose_file(n_1, l_1, n_2, l_2)
     mosh_bracket = 0
     if(abs(l_1 - l_2) <= abs(l - L)):
          x = abs(l - L)
     else:
          x = abs(l_1 - l_2)
     if(l_1 + l_2 <= l + L):
          y = l_1 + l_2
     else:
          y = l + L
     for lambd in xrange(x, y + 1):
          for i in xrange(0,len(z)):
               if(n == z[i][0] and l == z[i][1] and N == z[i][2] and L == z[i][3] and lambd == z[i][4]):
                    if(abs(z[i][6]) > 0.00000001):
                         mosh_bracket = z[i][6]
                         result = result + mosh_bracket*clebsch_gordan[lambd][l_1][l_2][m_1+6][m_2+6]*clebsch_gordan[lambd][l][L][m_l+6][M_L+6]
                         #conventie clebsch_gordan arrays zie matrices.py
     return result

#coefficient see thesis Maarten P74 eq 14.                            
def coefficient(S_,T, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N, L, M_L):           
     return (1./math.sqrt(2))*(1-(-1)**(S_+T+l))*clebsch_gordan_spin[S_][s_1][s_2]*clebsch_gordan_spin[T][t_1][t_2]*CG_MB_CG(n_1, l_1, m_1, n_2, l_2, m_2, n,l , m_l, N, L, M_L)
                    

def total_twobody(P, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, i):
     result = 0
     for N in xrange(0, 2*n_1 + l_1 +2*n_2 + l_2 + 1):
          for L in xrange(0, l_1 + l_2 + 1):
               a = Rad_wavefunc(P, v_mom(M_average,i), N, L)*norm(v_mom(M_average,i), N, L)
               for N_ in xrange(0, 2*n_1 + l_1 +2*n_2 + l_2 + 1):
                    b = Rad_wavefunc(P, v_mom(M_average,i), N_, L)*norm(v_mom(M_average,i), N_, L)
                    for n in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
                         for l in xrange(0, l_1 + l_2 + 1):
                              if(2*n_1 + l_1 + 2*n_2 + l_2  == 2*n + l + 2*N + L and 2*n_1 + l_1 + 2*n_2 + l_2 == 2*n + l + 2*N_ + L):
                                   for m_l in xrange(-l, l+1):
                                        for M_L in xrange(-L,L+1):
                                             for S_ in xrange(0,2):
                                                  for T in xrange(0,2):
                                                       result = result + coefficient(S_, T, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N, L, M_L)*coefficient(S_,T, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N_, L, M_L)*a*b
     return (1.0/(A[i]*(A[i]-1)))*result

shell = [0,1,2]
orbital = ["s","p","d","f"]
 
#Opvulling  van de schillen van Harmonische Oscillator   
def opvulling(P,i):
     result = 0
     for N_1 in xrange(0,shell[i]+1):
          if(N_1 % 2 == 0):
               start_1 = 0
          else:
               start_1 = 1
          for l_1 in xrange(start_1, N_1 + 1,2):
               n_1 = int((N_1 - l_1)/2.0)
               for m_1 in xrange(-l_1, l_1 +1):
                    for N_2 in xrange(0,shell[i]+1):
                         if(N_2 % 2 == 0):
                              start_2 = 0
                         else:
                              start_2 = 1
                         for l_2 in xrange(start_2, N_2 + 1, 2):
                              n_2 = int((N_2 - l_2)/2.0)
                              for m_2 in xrange(-l_2, l_2 + 1):
                                   for s_1 in xrange(0,2):
                                        for s_2 in xrange(0,2):
                                             for t_1 in xrange(0,2):
                                                  for t_2 in xrange(0,2):
                                                       if(n_1 != n_2 or l_1 != l_2 or m_1 != m_2 or s_1 != s_2 or t_1 != t_2):
                                                            result = result + total_twobody(P,n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, i)
     return result
     

def format_(value):
    return "%.8f" % value


step = 0.1
upperlimit = 25
for j in xrange(0,2):
     f = open("{:s}.txt".format(nuc[j]), "w")
     f.write("# mass number (A) = " + str(A[j]) +  '\n')
     f.write("# number of protons (Z) = " + str(Z[j]) +'\n')
     f.write("# k (1/fm)" + '\t' + "total two_body (fm^3)"  +  '\n')
     i = 0
     distr_array = []
     while (i < upperlimit):
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

