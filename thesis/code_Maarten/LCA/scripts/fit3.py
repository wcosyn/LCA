#!/usr/bin/python
# Gaussian fit to Maxwell - Boltzman distribution in outputfile (pdf) 
# Use python fit.py file outputfile ycolumn power_min_y_as power_max_y_as
import sys
import numpy as np
import pylab as pl
import math as m
from scipy import optimize


data = np.loadtxt(sys.argv[1])

fig = pl.figure(figsize=(8,6))


Xp = data[1:,0] 
#print Xp 
Xm = -1.*Xp
Xr = Xm[::-1].copy()
s = Xr.size
X = np.zeros( 2*s-1 )
X[:s-1] = Xr[:s-1]
X[s-1:] = Xp[0:s]
#X = np.sqrt(2)*X
#print X

Yf = data[1:,sys.argv[2]]
Yp = Yf
#print Yp
Ym=Yp
Yr = Ym[::-1].copy()
s = Yr.size
Y = np.zeros( 2*s-1 )
Y[:s-1] = Yr[:s-1]
Y[s-1:] = Yp[0:s]

#print Y
x = sum(X*Y)/sum(Y)
width = m.sqrt(abs(sum((X-x)**2*Y)/sum(Y)))

max = float(sys.argv[3])

fit = lambda x, width : max* np.sqrt(2/np.pi) / width**3 * np.exp( -(x)**2 / (2*width**2) )

#print Y

pl.xlabel( r' $P$ [$GeV$] ', fontsize='xx-large' )
pl.yscale('log')
p= pl.plot(Xp,fit(Xp, width), 'k-',lw=3)
pl.legend( [p], ['$\sigma_{CM}= %1.3f $'% width ] )
pl.plot(Xp, Yf, 'o', ms=5)
#pl.ylim((int(sys.argv[4]), int(sys.argv[5])))
#xlim((-2,2))


#print "succes? ", succes
print "results is ", width


#pl.savefig( sys.argv[2], format='pdf' )
pl.show()
