#cm two-body momentum for full shell nucleus
import math
from read import moshinsky
import scipy as sp
import numpy as np
from scipy import special
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_6j, wigner_9j
from scipy.integrate import quad
from collections import defaultdict


nuclide = raw_input("Nuclide? (He, Be, C, O, Al, Ca40, Ca48, Fe, Ag, Ar, Pb): ")
upperlimit = int(input("upperlimit momentum P (1/fm): "))
step = float(input("delta_P between data points(1/fm): "))

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
               break
     return result

def clebsch_gordan_spin(ms_1, ms_2, S):
     result = 0
     a = np.loadtxt("./cg/clebsch_spin.dat")
     for u in xrange(0,len(a)):
          if(ms_1 == a[u][0] and ms_2 == a[u][1] and S == a[u][2]):
               result = a[u][3]
               break
     return result
               

def CG_MB_CG(n_1, l_1, m_1, n_2, l_2, m_2, n,l , m_l, N, L, M_L):
     result = 0
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
                        

     
def hulp(n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, N, L, N_):
     key = (n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, N, L, N_)
     result = 0
     for n in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
          for l in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
               if(2*n_1 + l_1 + 2*n_2 + l_2  == 2*n + l + 2*N + L and 2*n_1 + l_1 + 2*n_2 + l_2 == 2*n + l + 2*N_ + L):
                    for m_l in xrange(-l, l+1):
                         for M_L in xrange(-L,L+1):
                              for S_ in xrange(0,2):
                                   for T in xrange(0,2):
                                        result = result + coefficient(S_,T, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N, L, M_L)*coefficient(S_ ,T, n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, n,l , m_l, N_, L, M_L)
     if(result != 0):
          cc[key] = result
     return result
   
     
def twobody(P,n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, i):
     result = 0
     for N in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
          for L in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
               a = Rad_wavefunc(P, v_mom(2*M_average, i), N, L)*norm(v_mom(2*M_average, i), N, L)
               for N_ in xrange(0, 2*n_1 + l_1 + 2*n_2 + l_2 + 1):
                    b = Rad_wavefunc(P, v_mom(2*M_average, i), N_, L)*norm(2*v_mom(M_average, i), N_, L)
                    if(P == 0):
                         result = result + (hulp(n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, N, L, N_)*a*b)
                    else:
                         key = (n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, N, L, N_)
                         result = result + a*b*cc[key]
     return (1.0/(A[i]*(A[i]-1)))*result    

def level(Q):
     j = 0
     s = 0
     while(s < Q):
          s = s + degen[j]
          j = j + 1
     return j


def opvulling(P,i):
     result = 0
     for t_1 in xrange(0,2):
          if(t_1 == 0):
               Q_1 = N[i]
          elif(t_1 == 1):
               Q_1 = Z[i]
          for lvl_1 in xrange(0, level(Q_1)):
               n_1 = radial[lvl_1]
               l_1 = orbital[lvl_1]
               a = float(degen[lvl_1])/ float(2*(2*l_1 + 1))
               if(nos[lvl_1]-Q_1 > 0):
                    b = float(Q_1-nos[lvl_1-1])/float(degen[lvl_1])
               else:
                    b = 1
               for t_2 in xrange(0,2):
                    if(t_2 == 0):
                         Q_2 = N[i]
                    elif(t_2 == 1):
                         Q_2 = Z[i]
                    for lvl_2 in xrange(0,level(Q_2)):
                         n_2 = radial[lvl_2]
                         l_2 = orbital[lvl_2]
                         if(t_1 == t_2 and n_1 == n_2 and l_1 == l_2 and degen[lvl_2] == degen[lvl_1]):
                              c = float(degen[lvl_2]-1)/ ((2*(2*l_2 + 1))-1)
                         else:
                              c = float(degen[lvl_2])/(2*(2*l_2 + 1))
                         if(nos[lvl_2]-Q_2 > 0):
                              if(t_1 == t_2 and n_1 == n_2 and l_1 == l_2 and degen[lvl_2] == degen[lvl_1]):
                                   d = float(Q_2-nos[lvl_2-1]-1)/float(degen[lvl_2]-1)
                              else:
                                   d = float(Q_2-nos[lvl_2-1])/float(degen[lvl_2])
                         else:
                              d = 1   
                         for m_1 in xrange(-l_1, l_1 + 1):
                              for m_2 in xrange(-l_2, l_2 + 1):
                                   for s_1 in xrange(0,2):
                                        for s_2 in xrange(0,2):
                                             if(t_1 != t_2 or n_1 != n_2 or l_1 != l_2 or m_1 != m_2 or s_1 != s_2 or degen[lvl_2] != degen[lvl_1]):
                                                  result = result + a*b*c*d*twobody(P,n_1, l_1, m_1, s_1, t_1, n_2, l_2, m_2, s_2, t_2, i)
     return result       

def format_(value):
    return "%.8f" % value


upperlimit = upperlimit/step


f = open("./cm/{:s}_tb_cm.txt".format(nuc[NUC]), "w")
f.write("# mass number (A) = " + str(A[NUC]) +  '\n')
f.write("# number of protons (Z) = " + str(Z[NUC]) +'\n')
f.write("# P (GeV)" + '\t' + "cm two_body (GeV^-3)"  +  '\n')
j = 0
distr_array = []
cc= defaultdict(int)
while (j <= upperlimit):
	print j
	c = opvulling(step*j, NUC)
	distr_array.append(c*(step*j)**2)
	txt = str(step*j) + '\t' + str(c) 
	f.write(txt)
	f.write('\n')
	j = j+1
nor = sp.integrate.trapz(distr_array,x=None,dx = step)
f.write("#normalisation = " + str(nor) +'\n') 
f.close()

