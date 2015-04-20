#calculate average kinetic energy of a nucleon
import math
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad
from collections import defaultdict


nuc = ["He","Be","C","O","Al","Ca40","Ca48","Fe","Ag","Ar","Pb"]
A   = [4,9,12,16,27,40,48,56,108,40,208]
Z   = [2,4,6,8,13,20,20,26,47,18,82]
N = [m - n for m,n in zip(A,Z)] #number of neutrons
     
hbarc = 197.327 #(MeV*fm)
M_n = 939.565 #mass neutron (MeV) [Particle Data Group]
M_p = 938.272 #mass proton (MeV) [Particle Data Group]
M_average = (M_n + M_p)/2.

f = open("kinetic_energy.txt", "w")
f.write('\t' + "p" + '\t'+ '\t'+ '\t' + "n" + '\n')
for i in nuc:
     fn = "./data/" + i + "_ob_mf.txt"
     data = np.loadtxt(fn)
     k = data[:,0]
     d_n = data[:,1]
     d_p = data[:,2]
     nor_n = sp.integrate.trapz(k*k*d_n, x=k)
     nor_p = sp.integrate.trapz(k*k*d_p, x=k)
     k_n = sp.integrate.trapz(k*k*k*k*d_n, x=k)/nor_n
     k_p = sp.integrate.trapz(k*k*k*k*d_p, x=k)/nor_p
     kin_n = (hbarc**2)*k_n/(2*M_n)
     kin_p = (hbarc**2)*k_p/(2*M_p)
     f.write(i + " : " + str(kin_p)+ '\t' + str(kin_n) +'\n')
f.close()
     
     
