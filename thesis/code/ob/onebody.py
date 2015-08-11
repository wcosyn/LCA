#one-body momentum 
import math
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad
from obStates import vulOp
import ctypes

#import radial wave function in c-code--------------------------------------------------------------------------------------------------------------------------------------------

integral = ctypes.cdll.LoadLibrary("../hulp/hoWave.so")
integral.nu.restype = ctypes.c_double
integral.nu.argtype =  ctypes.c_int
integral.radialWave.restype = ctypes.c_double
integral.radialWave.argtypes =  [ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_double]

#--------------------------------------------------------------------------------------------------------------------------------------


upperlimit = int(input("upperlimit momentum k (fm): "))
step = float(input("delta_k between data points(fm): "))

nuc = ["He","Be","C","O","Al","Ar","Ca40","Ca48","Fe","Ag","Pb"]
A   = [4,9,12,16,27,40,40,48,56,108,208]
Z   = [2,4,6,8,13,18,20,20,26,47,82]
N = [m - n for m,n in zip(A,Z)] #number of neutrons

Omega = []
for i in range(len(A)):
     Omega.append((45.0*(A[i]**(-1./3.)))-(25.0*(A[i]**(-2./3.))))
     
hbarc = 197.327 #(MeV*fm)
M_n = 939.565 #mass neutron (MeV) [Particle Data Group]
M_p = 938.272 #mass proton (MeV) [Particle Data Group]
M_average = (M_n + M_p)/2.
M_n = M_average
M_p = M_average

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

def max_l(i): 
     qq = level(N[i])
     nrgy = (2*radial + orbital)
     return nrgy[0:qq].max() + 1


def opvulling(k,t,i):
     result = 0
     ntrn = vulOp(i,0)
     prtn = vulOp(i,1)
     if(t==0):
          fill_1 = ntrn
          M = M_n
     else:     
          fill_1 = prtn
          M = M_p
     for N_1 in xrange(0,max_l(i)+1):
          if(N_1 % 2 == 0):
               n_1max = N_1+1
          else:
               n_1max = N_1
          for n in xrange(0,n_1max):
               l = N_1 - 2*(n)
               correctie_1 = fill_1[n][l]
               if(l>= 0):
                    result = result + correctie_1*(integral.radialWave(k,n,l,integral.nu(A[i])))**2
     return result       

def format_(value):
    return "%.3f" % value

def sigma(m1,m2):
     return np.sqrt(m2-m1**2)

def skewness(m1,m2,m3,sigma):
     return (m3-3*m1*m2+2*(m1**3))/(sigma**3)
     
def kurtosis(m1,m2,m3,m4,sigma):
     return (m4-4*m1*m3+6*m1**2*m2-3*m1**4)/(sigma**4) - 3.0 

file = open("ob_properties.txt", 'w')
file.write("#nucl" + '\t' + "sig2" + '\t' + "kurt" + '\t' + "skew" + '\n')
upperlimit = upperlimit/step

for i in xrange(0,2):
     f = open("./{:s}_ob.txt".format(nuc[i]), "w")
     f.write("# mass number (A) = " + str(A[i]) +  '\n')
     f.write("# number of protons (Z) = " + str(Z[i]) +'\n')
     f.write("#units k = 1/fm , units n_1(k) = fm^3"+'\n')
     f.write("# k (1/fm)" + '\t' + "n" +  '\t' + "p" +  '\t' + "tot"  + '\n')
     j = 0
     distr_array = []
     first_moment = []
     second_moment = []
     third_moment = []
     fourth_moment = []
     while (j <= upperlimit):
          print j
          P = opvulling(step*j, 1, i)
          Neu = opvulling(step*j, 0, i)
          tot = (Neu + P)/float(A[i])
          distr_array.append(tot*(step*j)**2)
          first_moment.append(tot*(step*j)**3)
          second_moment.append(tot*(step*j)**4)
          third_moment.append(tot*(step*j)**5)
          fourth_moment.append(tot*(step*j)**6)
          txt = str(step*j) + '\t' + str(Neu) + '\t' + str(P) + '\t' + str(tot)
          f.write(txt)
          f.write('\n')
          j = j+1
     nor = sp.integrate.trapz(distr_array,x=None,dx = step)
     one = sp.integrate.trapz(first_moment,x=None,dx = step)
     two = sp.integrate.trapz(second_moment,x=None,dx = step)
     three = sp.integrate.trapz(third_moment,x=None,dx = step)
     four = sp.integrate.trapz(fourth_moment,x=None,dx = step)
     f.write("#normalisation = " + str(nor) +'\n') 
     f.write("#eerste moment n(k)k**2 = " + str(one) +'\n') 
     f.write("#2e moment n(k)k**2 = " + str(two) +'\n') 
     f.write("#3e moment n(k)k**2 = " + str(three) +'\n') 
     f.write("#4e moment n(k)k**2 = " + str(four) +'\n')
     sig = sigma(one,two)
     f.write("#sigma**2 n(k)k**2 = " + str(sig**2) +'\n')
     f.write("#kurtosis n(k)k**2 = " + str(kurtosis(one,two,three,four,sig)) +'\n')
     f.write("#skewness n(k)k**2 = " + str(skewness(one,two,three,sig)) +'\n') 
     file.write("$nuclide[" + str(A[i]) + "][]{" + nuc[i] + "}$" + '\t' + "$" + format_(sig*hbarc)  +  "$" + '\t' +  "$" + format_(kurtosis(one,two,three,four,sig)) +  "$" + '\t' +  "$" + format_(skewness(one,two,three,sig))+  "$" +'\n')
     f.close()

