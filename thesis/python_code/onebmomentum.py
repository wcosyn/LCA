#calculate one-body momentum distribution

import math 
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


#returns HO constant \nu for different mass
def v(M, index):
     return M*Omega[index] / (hbarc**2)

def v_mom(M, index):
     return (hbarc**2) / (Omega[index]*M) 

#returns unnormalized radial wavefunction calculated with constant a 
def Rad_wavefunc(k, a, n, l):
	t = a*(k**2)
	return (k**l)*sp.special.eval_genlaguerre(n, l+0.5, t)*np.exp(-(a/2.)*k**2)

def norm(a, n, l):
	return np.sqrt((2*math.factorial(n)*a**(l+(3./2.)))/math.gamma(n+l+(3./2.)))
     


degen   = [2,4,2,6,2,4,8,4,6,2,10,8,6,4,2,12,10,8,6,4,2,14,10,6,12,8,2,4]
orbital = [0,1,1,2,0,2,3,1,3,1,4,4,2,2,0,5,5,3,3,1,1,6,4,2,6,4,0,2]
radialQ = [0,0,0,0,1,0,0,1,0,1,0,0,1,1,2,0,0,1,1,2,2,0,1,2,0,1,3,2] 

def onebodym(k, Q, M, i):
     j = 0
     s = 0
     result = 0
     level = 0
     f = 0
     while(s < Q):
          s = s + degen[j]
          j = j + 1
     counter = 0
     while(level <= j):
          if ((counter + degen[level]) <= Q):
               f = degen[level]
          else:
               f = Q - counter
          result = result + f*(norm(v_mom(M, i), radialQ[level], orbital[level])*Rad_wavefunc(k, v_mom(M, i), radialQ[level], orbital[level]))**2
          level = level + 1
          counter = counter + f
     return result
     
          
          
          
     
     
	
def one_b_mom_dist(k, Q, M, i):
     f = 0	
     n = 0
     counter = 0
     result = 0
     l = 0
     d = 0
     N = 0 
     while (counter < Q):
          l = N
          while (l >= 0 and counter < Q):
                    d = 2*((2*l)+1)
                    n = (N-l)/2.0
                    if ((counter + d) <= Q):
                         f = 1
                    else:
                         f = float(Q - counter)/float(d)
                    result = result + d*f*(norm(v_mom(M, i), n, l)*Rad_wavefunc(k, v_mom(M, i), n, l))**2
                    counter = counter + int(d*f)
                    l = l - 2
          N = N + 1
     return result

                    
                            
def distribution(k, a, d, b, i ):
     return (onebodym(k, d, M_n, i ) + onebodym(k, b, M_p, i))*a/a

def integrand(k, a, d, b, i):
        return distribution(k,a, d, b, i)*(k**2)



for j in range(len(nuc)):
     f = open("{:s}_mfA.txt".format(nuc[j]), "w")
     f.write("# mass number (A) = " + str(A[j]) +  '\n')
     f.write("# number of protons (Z) = " + str(Z[j]) +'\n')
     f.write("# k (1/fm)" + '\t' + "mf (fm^3)" + '\n')
     i = 0
     while (i < 25):
          txt = str(0.1*i) + '\t' + str(distribution(0.1*i, A[j], N[j], Z[j], j))
          f.write(txt)
          f.write('\n')
          i = i+1
     nor, err = sp.integrate.quad(integrand, 0, 10, args = (A[j], N[j], Z[j], j))
     f.write("#normalisation = " + str(nor) + '\n') 
     f.write("#integration error = " + str(err))
     f.close()

