#!/usr/bin/python
import sys
import numpy as np
import pylab as pl
from scipy import special
import math as m
import matplotlib.cm as cm
from matplotlib import rc
rc('text', usetex=True)
rc('font', size=24)
#rc('font',family='serif')


def norm( nu, n, l):
  return m.sqrt( 2.* m.gamma(n+1)/ m.gamma(n+l+1.5) * m.pow( nu, l+1.5))

def ho( nu, n, l, r):
  return norm( nu, n, l)* np.power(r, l)* np.exp( -0.5*nu*r*r) * special.hyp1f1(-n, l+1.5, nu*r*r)* special.binom(n+l+0.5, n)


fig=pl.figure()
fig.subplots_adjust(bottom=0.15, left=0.15)

A= int(sys.argv[1])
n= int(sys.argv[2])
l= int(sys.argv[3])
j= int(sys.argv[4])
name= sys.argv[5]

hbaromega= 45* m.pow( A, -1./3.) - 25* m.pow( A, -2./3.)
nu= 938* hbaromega/ 197.327/ 197.327

X,Y = np.loadtxt("WF{:s}{:d}{:d}{:d}1.dat2".format(name, n, l, j), usecols=[0,3], unpack=True )

nHO, coef= np.loadtxt("wsexp{:d}{:d}{:d}{:d}1.dat2".format(A, n, l, j), unpack=True)

print nHO, coef


result= np.zeros_like(X)


sum_a= 0
for i in range( 0, len(nHO)):
  print nHO[i], coef[i]
  result+= coef[i]* ho( nu, nHO[i], l, X)
  sum_a+= coef[i]*coef[i]

print sum_a
print sum( coef*coef)





pl.plot( X, Y, '.', label=r'$\bm{^{12}}$ C', lw=3)
pl.plot( X, result, label="HO expansion", lw=3)
pl.xlabel("r [fm]")
pl.ylabel(r"$R^{WS}_{00\frac{1}{2}}(r) [$fm$^{-3/2}]$")
pl.legend(loc="upper right")



pl.show()




