#!/usr/bin/python
# Fit maxwel-boltzman distribution to the data from inputfile and plot results to pdf-file outputfile
# Use fit.py inputfile output-pdf-file ycolumn min_y_as max_y_as
import sys
import numpy as np
import pylab as pl
from scipy import optimize


data = np.loadtxt(sys.argv[1])

fig = pl.figure(figsize=(8,6))


X = data[:,0]
#print Xp 
#X = np.sqrt(2)*X
#print X

Y = data[:,sys.argv[3]]
Y = X*X*Y
#print Yp 
#print Y

fitfunc = lambda x, a, b: a* np.sqrt(2/np.pi) *x*x / b**3 * np.exp( -(x)**2 / (2*b**2) ) # Target function
rfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
p0 = [100, 0.170] # Initial guess for the parameters
#p1, succesm  = optimize.leastsq(errfunc, p0[:], args=(X, Y))

popt, pcov = optimize.curve_fit(fitfunc, X, Y, p0)
print 'best fit parameters ',popt
print 'cov matrix', pcov

pl.xlabel( "$ P [fm^{-1}]$" )
p= pl.plot(X,fitfunc(X, popt[0], popt[1]), 'k-',lw=3)
pl.legend( p, ["\sigma_{CM}= %1.3f \; fm^{-1} "%popt[1] ] )
pl.plot(X, Y, 'o', ms=5)
pl.ylim((int(sys.argv[4]), int(sys.argv[5])))
pl.yscale('log')




#pl.savefig( sys.argv[2], format='pdf' )
pl.show()
