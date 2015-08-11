#!/usr/bin/python
# GOOD FITTING FUNCTION USING LEASTSQR METHOD
# ARGS: FILE COLX COLY COLYERR
import sys
import io
import math as m
import numpy as np
import pylab as pl
from scipy import optimize

# Load Data
if len(sys.argv ) == 5 :
#  data = np.loadtxt(sys.argv[1], usecols= (int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])) )
  data = np.loadtxt(sys.argv[1], usecols= (int(sys.argv[2]), int(sys.argv[3]), 0, 5, 6, 7, 8 ) )

  # Fit Gaussian to function
  X = (data[0:,3]+ data[:,4]+ data[:,5]+ data[:,6])/data[:,2]*2

#  X = data[:,0]
  Y = data[:,0]
  ERRY = data[:,1]
   
  print X
  print Y
  print ERRY

  fitfunc = lambda x, p0, p1: p0 + x*p1  # Target function

  pinit = [1, 1] # Initial guess for the parameters

  pbest = optimize.curve_fit(fitfunc,X, Y, pinit, ERRY)

  bestparams = pbest[0]
  cov_x = pbest[1]

  print bestparams
  print cov_x
  print np.sqrt(cov_x)

  pl.title( sys.argv[1] )
  p= pl.plot(X,fitfunc(X, bestparams[0], bestparams[1] ), 'k-',lw=3)
  pl.errorbar(X, Y, yerr=ERRY )

  pl.show()
else :
  data = np.loadtxt(sys.argv[1], usecols= (int(sys.argv[2]), int(sys.argv[3])) )

  # Fit Gaussian to function
  X = np.log10(data[:,0])
  Y = np.log10(data[:,1])
  print X
  print Y

  fitfunc = lambda x, p0, p1: p0 + x*p1  # Target function

  pinit = [1, 1] # Initial guess for the parameters

  pbest = optimize.curve_fit(fitfunc,X, Y, pinit)

  bestparams = pbest[0]
  cov_x = pbest[1]

  print "bestfit", bestparams
  print cov_x
  print np.sqrt(cov_x)

  pl.title( sys.argv[1] )
#  pl.yscale('log')
#  pl.xscale('log')
  #p= pl.semilogx(X,fitfunc(X, bestparams[0], bestparams[1] ), 'k-',lw=3)
  p= pl.plot(X,fitfunc(X, bestparams[0], bestparams[1] ), 'k-',lw=3)
  pl.plot(X, Y )

  pl.show()


# Close the files


