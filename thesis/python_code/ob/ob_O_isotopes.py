#total two-body momentum for full shell nucleus
import math
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad

upperlimit = int(input("upperlimit momentum k (fm): "))
step = float(input("delta_k between data points(fm): "))

nuc = ["O12","O13","O14","O15","O16","O17","O18","O19","O20","O21","O22","O23","O24"]
A   = [12,13,14,15,16,17,18,19,20,21,22,23,24]
Z   = s[8,8,8,8,8,8,8,8,8,8,8,8,8]
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
     
def level(Q):
     j = 0
     s = 0
     while(s < Q):
          s = s + degen[j]
          j = j + 1
     return j


def opvulling(k,i):
     result = 0
     for t in xrange(0,2):
          if(t == 0):
               Q = N[i]
               M = M_average
          elif(t == 1):
               Q = Z[i]
               M = M_average
          for lvl in xrange(0, level(Q)):
               n = radial[lvl]
               l = orbital[lvl]
               a = float(degen[lvl])/ float(2*(2*l + 1))
               if(nos[lvl]-Q > 0):
                    b = float(Q-nos[lvl-1])/float(degen[lvl])
               else:
                    b = 1
               result = result + a*b*(2*(2*l + 1))*(norm(v_mom(M, i), n, l)*Rad_wavefunc(k, v_mom(M, i), n, l))**2
     return (1.0/A[i])*result       

def format_(value):
    return "%.8f" % value


upperlimit = upperlimit/step

for i in xrange(0,len(nuc)):
     f = open("./data/{:s}_ob_mf_Mav.txt".format(nuc[i]), "w")
     f.write("# mass number (A) = " + str(A[i]) +  '\n')
     f.write("# number of protons (Z) = " + str(Z[i]) +'\n')
     f.write("# number of neutrons (N) = " + str(N[i]) +'\n')
     f.write("# k (1/fm)" + '\t' + "n_1(k) (fm^3)"  +  '\n')
     j = 0
     distr_array = []
     while (j <= upperlimit):
          print j
          c = opvulling(step*j, i)
          distr_array.append(c*(step*j)**2)
          txt = str(step*j) + '\t' + str(c) 
          f.write(txt)
          f.write('\n')
          j = j+1
     nor = sp.integrate.trapz(distr_array,x=None,dx = step)
     f.write("#normalisation = " + str(nor) +'\n') 
     f.close()
