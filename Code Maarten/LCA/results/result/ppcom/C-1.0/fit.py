#!/usr/bin/python
# ARG: file	xcol	ycol	

from fitting import *

import sys
from scipy import optimize
from numpy import * 
import pylab as pl



# Load Data
data = loadtxt(sys.argv[1])
X= data[:,sys.argv[2]]
Y= data[:,sys.argv[3]]
dX= X[5]-X[4]

i= 0
for p in X :
  if p < 0.600: i+=1
  else: break

X= X[:i]
Y= Y[:i]

#Check integration 
integr= sum( X*X*Y*dX )
print "Integral is ", integr


# Begin of fit code.
# giving initial parameters
mu = 0.01
sigma = 0.1
height = 10
p= [sigma,height]

# define your function:
def f(x, s, h ): return h * exp(-0.5*((x)/s)**2)

# fit! (given that data is an array with the data to fit)
# do fit using Levenberg-Marquardt

popt, punc, rchi2, dof = general_fit( f, X, Y, p )


for i in range( len(p) ):
  print "%2i %12f +/- %10f"%(i,popt[i],punc[i])


print popt[0] * sqrt(2.), "+/-", punc[0] * sqrt(2)
# Plot

# plot data
# errorbar(freq, vr/v0, yerr=uvr, fmt='r+', label="Data")
# plot fit
pl.plot(X, Y, 'g-', label="Start")
pl.plot(X, f(X, popt[0], popt[1]), 'b-', label="Fit")
pl.yscale("log")
#xlabel("Frequency [Hz]")
#ylabel("Relative amplitude")
#legend(("data","fit"))
#legend()
pl.show()


