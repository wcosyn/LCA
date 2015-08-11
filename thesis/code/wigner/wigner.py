#Wignerdistributie voor algemene HO golffunctie
#imports---------------------------------------------------------------------------------------------------------------------------------------

import math
import scipy as sp
import numpy as np
import cmath
from scipy import special
from scipy import integrate
import sys
import ctypes

#import radial wave function in c-code--------------------------------------------------------------------------------------------------------------------------------------------

integral = ctypes.cdll.LoadLibrary("./c_wigner/integral.so")
integral.nu.restype = ctypes.c_double
integral.nu.argtype =  ctypes.c_int
integral.wigner_integrand.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_int]
integral.wigner_integrand.restype = ctypes.c_double
integral.theta.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double]
integral.theta.restype = ctypes.c_double

#inputs---------------------------------------------------------------------------------------------------------------------------------------

nuclide = sys.argv[1]
upperlimit = float(sys.argv[3])
lowerlimit = float(sys.argv[2])
step = float(sys.argv[4])
NN = int(sys.argv[5])
LL = int(sys.argv[6])
MM = int(sys.argv[7])


#nuclear parameters----------------------------------------------------------------------------------------------------------------------------
     
nuc = ["He","Be","C","O","Al","Ca40","Ca48","Fe","Ag","Ar","Pb"]
A   = [4,9,12,16,27,40,48,56,108,40,208]
Z   = [2,4,6,8,13,20,20,26,47,18,82]
N = [m - n for m,n in zip(A,Z)] #number of neutrons

NUC = 0
if nuclide in nuc:
	for i in xrange(0,len(nuc)):
		if(nuclide == nuc[i]):
			NUC = i
			break
               
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
     
#main----------------------------------------------------------------------------------------------------------------------------------------------

#norm wigner function
def norm_wigner(l,m):
     if(math.fabs(m)<=l):
          return ((-1.0)**(m))*(2*l+1)*math.factorial(l-m)/(2*(np.pi)*math.factorial(l+m))

def integrand(ksi,x_,n,l,m,r,k,i):
     return integral.wigner_integrand(ksi,x_,n,l,m,r,k,i);

def r_min(x):
     return 0

def r_max(x):
     return 300.0

def wigner(n,l,m,r,k,i):
     return integrate.dblquad(integrand,-1,1, r_min, r_max, args=(n,l,m,r,k,A[i]))[0]


      
#calculate on grid----------------------------------------------------------------------------------------------------------------------------------
    
interval = "_" + str(lowerlimit) + "_" + str(upperlimit)
upperlimit = int(upperlimit/step)
lowerlimit = int(lowerlimit/step)


f = open("./{:s}/{:s}_wigner_".format(nuc[NUC],nuc[NUC]) + str(NN)  + str(LL)+ str(MM) + "_" + interval + ".txt", "w")
ff = open("./{:s}/{:s}_ob_test_".format(nuc[NUC],nuc[NUC])   + str(NN)  + str(LL) + str(MM) + "_" + interval + ".txt", "w")
j = lowerlimit
NORM = norm_wigner(LL,MM)
while (j < upperlimit):
     K = j*step
     f.write('\n')
     print "k = " + str(K)
     w = 0.0
     distr = []
     result = 0
     while(w <= 6.0/step):
          R = w*step
          c = NORM*wigner(NN,LL,MM, R, K, NUC)
          if(R == 3.0):
               print "Half way!"
          txt = str(K) + '\t' + str(R) + '\t' +str(c) 
          distr.append(c*R*R)
          f.write(txt)
          f.write('\n')
          w = w + 1
     ff.write( str(K) + '\t' + str(integrate.trapz(distr,x=None,dx=step)) + '\n')
     j = j+1
f.close()
ff.close()
"""

kk_ = 2.0
ksi = 2.0
xx_ = 0.4
step = 0.00001
f = open("test_integrand.txt",'w')
for i in xrange(-200000,200001):
     point = i*step
     f.write(str(point) + '\t' + str(integral.theta(point,xx_,kk_,ksi)) + '\n')
print "done"
f.close()"""






