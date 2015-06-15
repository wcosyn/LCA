#!/usr/bin/python
import sys
import io
import math as m
import numpy as np
import pylab as pl
from scipy import optimize

# Load Data
data = np.loadtxt(sys.argv[1])
X= data[:,sys.argv[2]]
Y= data[:,sys.argv[3]]*X*X

dX= X[5]-X[4]
integr= sum( Y*dX )

fitfunc = lambda x, p0, p1: p0*x*x/p1/p1/p1*np.exp( -x*x / (2*p1**2) ) # Target function
p0 = [10, 0.15] # Initial guess for the parameters
max = Y.max()
#p1, succesm  = optimize.leastsq(errfunc, p0[:], args=(X, Y))

pbest = optimize.curve_fit(fitfunc, X, Y, p0) 

bestparams = pbest[0]
cov_x = pbest[1]


#pl.title( sys.argv[1] )
pl.xlabel( r' $P$ [$MeV$] ', fontsize='xx-large' )
p= pl.plot(X,fitfunc(X, bestparams[0], bestparams[1] )/max, 'k-',lw=3)
pl.legend( [p], ['\sigma_{CM}= %(two)1.3f$'%{  "two" : bestparams[1]  } ] )
pl.plot(X, Y/max, 'o', ms=5)
pl.yscale('log')
#pl.ylim( 1e-3, 1e0 )

print 'best fit parameters ',bestparams
#print cov_x
print sys.argv[1], sys.argv[3], "Width", bestparams[1], "+-", m.sqrt(cov_x[1][1])
print "Total integral ", integr
print bestparams


pl.show()
# Close the files
