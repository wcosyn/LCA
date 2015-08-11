#total two-body momentum for full shell nucleus
import math
from read import moshinsky
import scipy as sp
import numpy as np
from scipy import special
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_6j, wigner_9j
from scipy.integrate import quad
from collections import defaultdict
import sys
import ctypes

ho_lib = ctypes.cdll.LoadLibrary("./HO_rad/ho_rad.so")
ho_lib.R.restype = ctypes.c_double
ho_lib.R.argtypes = [ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_double]
ho_lib.nu.restype = ctypes.c_double
ho_lib.nu.argtype =  ctypes.c_int


nuclide = sys.argv[1]
upperlimit = float(sys.argv[2])
step = float(sys.argv[3])

nuc = ["He","Be","C","O","Al","Ca40","Ca48","Fe","Ag","Ar","Pb"]
A   = [4,9,12,16,27,40,48,56,109,40,208]
Z   = [2,4,6,8,13,20,20,26,47,18,82]
N = [m - n for m,n in zip(A,Z)] #number of neutrons

NUC = 0
if nuclide in nuc:
	for i in xrange(0,len(nuc)):
		if(nuclide == nuc[i]):
			NUC = i
			break
else:
	print "nuclide not found"
Omega = []
for i in range(len(A)):
     Omega.append((45.0*(A[i]**(-1./3.)))-(25.0*(A[i]**(-2./3.))))
     
hbarc = 197.327 #(MeV*fm)
M_n = 938. #mass neutron (MeV) [Particle Data Group]
M_p = 938. #mass proton (MeV) [Particle Data Group]
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
     

#returns HO constant \nu for different mass
def v(M, index):
     return M*Omega[index] / (hbarc**2)  
     
#returns HO constant \nu in momentum space for different mass
def v_mom(M, index):
     return (hbarc**2) / (Omega[index]*M) 

#returns unnormalized radial wavefunction for HO calculated with constant a 
def Rad_wavefunc(k, a, n, l):
	t = a*(k**2)
	return (k**l)*special.eval_genlaguerre(n, l+0.5, t)*np.exp(-(a/2.)*(k**2))
     
#norm radial wave function
def norm(a, n, l):
	return np.sqrt((2*math.factorial(n)*(a**(l+3./2.)))/math.gamma(n+l+(3./2.)))

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
     
def clebsch_gordan(l_1, m_1, l_2, m_2 , L):
     result = 0
     filename = "./cg/cg" + str(l_1) + str(l_2) + str(L) + ".dat"
     a = np.loadtxt(filename)
     for u in xrange(0,len(a)):
          if(m_1 == a[u][0] and m_2 == a[u][1]):
               result = a[u][2]
     return result

def clebsch_gordan_spin(ms_1, ms_2, S):
     result = 0
     a = np.loadtxt("./cg/clebsch_spin.dat")
     for u in xrange(0,len(a)):
          if(ms_1 == a[u][0] and ms_2 == a[u][1] and S == a[u][2]):
               result = a[u][3]
     return result
               
def CG_MB_CG(n_1, l_1, m_1, n_2, l_2, m_2, n,l , m_l, N, L, M_L):
     key_ = (n_1, l_1, m_1, n_2, l_2, m_2, n,l,m_l, N, L, M_L)
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
     if(result !=0):
          dd[key_] = result 
     return result

def coefficient(S_, T, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N, L, M_L): 
     key_ =  (n_1, l_1, m_1, n_2, l_2, m_2, n,l,m_l, N, L, M_L)
     fact = 0   
     if(dd[key_] == 0 ):
        fact = CG_MB_CG(n_1, l_1, m_1, n_2, l_2, m_2, n,l , m_l, N, L, M_L)
     else:
        fact = dd[key_]
     return (1.0/math.sqrt(2))*(1.0-(-1.0)**(l + S_ + T))*(clebsch_gordan_spin(s_1, s_2, S_))*(clebsch_gordan_spin( t_1, t_2, T))*fact
                        

     
def hulp(n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,n_,l):
     key = (n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,n_,l)
     result = 0
     for N in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
          for L in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
               if(2*n_1 + l_1 + 2*n_2 + l_2  == 2*n + l + 2*N + L and 2*n_1 + l_1 + 2*n_2 + l_2 == 2*n_ + l + 2*N + L):
                    for m_l in xrange(-l, l+1):
                         for M_L in xrange(-L,L+1):
                              for S_ in xrange(0,2):
                                   for T in xrange(0,2):
                                        result = result + coefficient(S_,T, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N, L, M_L)*coefficient(S_,T, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n_,l , m_l, N, L, M_L)
     if(result != 0):
          cc[key] = result
     return result
   
     
def twobody(P,n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, i,l):
     result = 0
     for n in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
          a = ho_lib.R(P,n,l,ho_lib.nu(A[i]))
          for n_ in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
               b = ho_lib.R(P,n_,l,ho_lib.nu(A[i]))
               if(P == 0):
                    result = result + (hulp(n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,n_,l)*a*b)
               else:
                    key = (n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2,n,n_,l)
                    result = result + a*b*cc[key]
     return result    

def level(Q):
     j = 0
     s = 0
     while(s < Q):
          s = s + degen[j]
          j = j + 1
     return j


def vulOp(i,t):
     opv =  [[0 for x in range(5)] for x in range(9)] 
     result = 0
     if(t == 0):
          Q = N[i]
     elif(t == 1):
          Q = Z[i]
     for lvl in xrange(0, level(Q)):
          n = radial[lvl]
          l = orbital[lvl]
          if(Q < nos[lvl]):
               opv[n][l] = opv[n][l] + Q-nos[lvl-1]
          else:
               opv[n][l] = opv[n][l] + degen[lvl]
     return opv
                         

def opvulling(P,i,l):
     result = 0
     ntrn = vulOp(i,0)
     prtn = vulOp(i,1)
     for t_1 in xrange(0,2):
          if(t_1==0):
               fill_1 = ntrn
          else:     
               fill_1 = prtn 
          for N_1 in xrange(0,max_l(i)+1):
               if(N_1 % 2 == 0):
                    n_1max = N_1+1
               else:
                    n_1max = N_1
               for n_1 in xrange(0,n_1max):
                    l_1 = N_1 - 2*(n_1)
                    degen_1 = 2*(2*l_1+1) 
                    correctie_1 = float(fill_1[n_1][l_1])/float(degen_1)
                    for t_2 in xrange(0,2):
                         if(t_2==0):
                              fill_2 = ntrn
                         else:     
                              fill_2 = prtn 
                         for N_2 in xrange(0,max_l(i)+1):
                              if(N_2 % 2 == 0):
                                   n_2max = N_2+1
                              else:
                                   n_2max = N_2
                              for n_2 in xrange(0,n_2max):
                                   l_2 = N_2 - 2*(n_2)
                                   degen_2 = 2*(2*l_2+1)
                                   if(n_1 == n_2 and l_1 == l_2 and t_1 == t_2):
                                        correctie_2 = float(fill_2[n_2][l_2]-1)/float(degen_2-1)
                                   else:
                                        correctie_2 = float(fill_2[n_2][l_2])/float(degen_2)
                                   if(l_1 >= 0 and l_2 >= 0):
                                        for m_1 in xrange(-l_1, l_1+1):
                                             for s_1 in xrange(0,2):
                                                  for m_2 in xrange(-l_2,l_2 + 1):
                                                       for s_2 in xrange(0,2):
                                                            result = result + correctie_1*correctie_2*twobody(P,n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, i,l)
     return (1.0/(A[i]*(A[i]-1)))*result

def format_(value):
    return "%.8f" % value



def max_l(): 
     qq = level(N[NUC])
     nrgy = (2*radial + orbital)*2
     return nrgy[0:qq].max() + 1 

def sigma(m1,m2):
     return np.sqrt(m2-m1**2)

def skewness(m1,m2,m3,sigma):
     return (m3-3*m1*m2+2*(m1**3))/(sigma**3)
     
def kurtosis(m1,m2,m3,m4,sigma):
     return (m4-4*m1*m3+6*m1**2*m2-3*m1**4)/(sigma**4) - 3.0 


upperlimit = upperlimit/step
dd = defaultdict(float)
for l in xrange(0,2*max_l(NUC)+1):
     fn = "./relative/new_results/{:s}_tb_rel_".format(nuc[NUC]) + str(l)  + ".txt"
     f = open(fn, "w")
     f.write("# mass number (A) = " + str(A[NUC]) +  '\n')
     f.write("# number of protons (Z) = " + str(Z[NUC]) +'\n')
     f.write("#relative orbital quantum number l = " + str(l) + '\n')
     f.write("# k (1/fm)" + '\t' + "cm two_body (fm^3)"  +  '\n')
     j = 0
     distr_array = []
     first_moment = []
     second_moment = []
     third_moment = []
     fourth_moment = []
     cc= defaultdict(float)
     while (j <= upperlimit):
          mom = j*step
          print nuclide + '\t' + str(mom)
          c = opvulling(mom, NUC,l)
          distr_array.append(c*(mom)**2)
          first_moment.append(c*(mom)**3)
          second_moment.append(c*(mom)**4)
          third_moment.append(c*(mom)**5)
          fourth_moment.append(c*(mom)**6)
          txt = str(mom) + '\t' + str(c) 
          f.write(txt)
          f.write('\n')
          j = j+1
     nor = sp.integrate.trapz(distr_array,x=None,dx = step)
     one = sp.integrate.trapz(first_moment,x=None,dx=step)
     two = sp.integrate.trapz(second_moment,x=None,dx=step)
     three = sp.integrate.trapz(third_moment,x=None,dx=step)
     four = sp.integrate.trapz(fourth_moment,x=None,dx=step)
     f.write("#normalisation = " + str(nor) +'\n') 
     sig = sigma(one,two)
     f.write("#sigma**2 n(k)k**2 = " + str(sig**2) +'\n')
     f.write("#kurtosis n(k)k**2 = " + str(kurtosis(one,two,three,four,sig)) +'\n')
     f.write("#skewness n(k)k**2 = " + str(skewness(one,two,three,sig)) +'\n') 
     f.close()
