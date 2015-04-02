#total two-body momentum for full shell nucleus
import math
from read import moshinsky
import scipy as sp
import numpy as np
from scipy import special
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_6j, wigner_9j
from scipy.integrate import quad


nuc = ["He","Be","C","O","Al","Ca40","Ca48","Fe","Ag","Ar","Pb"]
A   = [4,9,12,16,27,40,48,56,109,40,208]
Z   = [2,4,6,8,13,20,20,26,47,18,82]
N = [m - n for m,n in zip(A,Z)] #number of neutrons

Omega = []
for i in range(len(A)):
     Omega.append((45.0*(A[i]**(-1./3.)))-(25.0*(A[i]**(-2./3.))))
     
hbarc = 197.327 #(MeV*fm)
M_n = 939.565 #mass neutron (MeV) [Particle Data Group]
M_p = 938.272 #mass proton (MeV) [Particle Data Group]
M_average = (M_n + M_p)/2.

degen   = np.array([2,4,2,6,2,4,8,4,6,2,10,8,6,4,2,12,10,8,6,4,2,14,10,6,12,8,2,4])
orbital = np.array([0,1,1,2,0,2,3,1,3,1,4,4,2,2,0,5,5,3,3,1,1,6,4,2,6,4,0,2])
radial = np.array([0,0,0,0,1,0,0,1,0,1,0,0,1,1,2,0,0,1,1,2,2,0,1,2,0,1,3,2])
nos = []
for i in xrange(0,len(degen)):
     s = 0
     for j in xrange(0,i+1):
           s = s + degen[j]
     nos.append(s)
     

def clebsch_gordan(l_1, m_1, l_2, m_2 , L):
     result = 0
     filename = "cg" + str(l_1) + str(l_2) + str(L) + ".dat"
     a = np.loadtxt(filename)
     for u in xrange(0,len(a)):
          if(m_1 == a[u][0] and m_2 == a[u][1]):
               result = a[u][2]
               break
     return result

def clebsch_gordan_spin(ms_1, ms_2, S):
     result = 0
     a = np.loadtxt("clebsch_spin.dat")
     for u in xrange(0,len(a)):
          if(ms_1 == a[u][0] and ms_2 == a[u][1] and S == a[u][2]):
               result = a[u][3]
               break
     return result

               

def CG_MB_CG(n_1, l_1, m_1, n_2, l_2, m_2, n,l , m_l, N, L, M_L):
     result = 0
     mosh_bracket = 0
     if(math.fabs(l_1 - l_2) <= math.fabs(l - L)):
          x = int(math.fabs(l - L))
     else:
          x = int(math.fabs(l_1 - l_2))
     if(l_1 + l_2 <= l + L):
          y = l_1 + l_2
     else:
          y = l + L
     if(m_1 + m_2 == m_l + M_L):
          for lambd in xrange(x, y + 1):
               result = result + (moshinsky(n_1, l_1, n_2, l_2, n, l, N, L, lambd)*(clebsch_gordan(l_1, m_1, l_2, m_2, lambd))*(clebsch_gordan(l, m_l, L, M_L, lambd)))
     return result

def coefficient(S_, T, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N, L, M_L):           
     return (1.0/math.sqrt(2))*(1.0-(-1.0)**(l + S_ + T))*CG_MB_CG(n_1, l_1, m_1, n_2, l_2, m_2, n,l , m_l, N, L, M_L)*(clebsch_gordan_spin(s_1, s_2, S_))*(clebsch_gordan_spin( t_1, t_2, T))


for S_ in xrange(0,2):
     for T in xrange(0,2):
          for n_1 in xrange(0,4):
               for l_1 in xrange(0,7):
                    for n_2 in xrange(0,4):
                         for l_2 in xrange(0,7):
                              for n in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
                                   for l in xrange(0, l_1 + l_2 + 1):
                                        for N in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
                                             for L in xrange(0,l_1 + l_2 + 1):
                                                  file_name = "coef" + str(1000000000*S_ + 100000000*T + 10000000*n_1 + 1000000*l_1 + 100000*n_2 + 10000*l_2 + 1000*n + 100*l + 10*N + 1*L) + ".dat"
                                                  f = open(file_name, "w")
                                                  for m_1 in xrange(-l_1, l_1 + 1):
                                                       for m_2 in xrange(-l_2, l_2 + 1):
                                                            for m_l in xrange(-l, l + 1):
                                                                 for M_L in xrange(-L, L + 1):
                                                                      for s_1 in xrange(0,2):
                                                                           for s_2 in xrange(0,2):
                                                                                for t_1 in xrange(0,2):
                                                                                     for t_2 in xrange(0,2):
                                                                                          txt = str(m_1) + '\t' + str(m_2) + '\t' + str(s_1) + '\t' + str(s_2) + '\t' + str(t_1) + '\t' + str(t_2) + '\t' +  str(coefficient(S_, T, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N, L, M_L))
                                                                                f.write(txt)
                                                                                f.write('\n')
                                                  f.close()


                                                  
                              
                                                  
                              