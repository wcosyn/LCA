#calcnorm maarten
import numpy as np
import scipy as sp
import math
from scipy import special
from scipy.integrate import quad

nuclide = int(input("Nuclide? (He, C, Al, Fe, Pb): "))

nuc = ["He","C","Al","Fe","Pb"]
fn = nuc[nuclide] + "_" + "cm"  + "_" + "maarten" + ".txt"
data = np.loadtxt(fn)
dist= data[:,1]
k = data[:,0]
kquad = np.multiply(k,k)
integrand = np.multiply(kquad,dist)
nor = sp.integrate.trapz(integrand,x=data[:,0])

data[:,1] = data[:,1]/nor

fn = nuc[nuclide] + "_cm" + "_maarten_normed.txt"
f = open(fn,"w")
f.write("#k"+ '\t' + "cm_distr" + '\n')
for i in xrange(0, len(data[:,1])):
     f.write(str(data[i,0]) + '\t' +  str(data[i,1]) + '\n')
f.close()

print nor

k2= data[:,0]
kquad2 = np.multiply(k2,k2)
dist2= data[:,1]
integrand2 = np.multiply(kquad2,dist2)
nor2 = sp.integrate.trapz(integrand2,x=data[:,0])
print nor2
     
