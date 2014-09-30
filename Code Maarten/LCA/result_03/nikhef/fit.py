#!/usr/bin/python
# ARG: file	xcol	ycol	 y_err

from fitting import *

import sys
from scipy import optimize
from numpy import * 
import pylab as pl



# Load Data
data = loadtxt(sys.argv[1])
X= data[:,sys.argv[2]]
Y= data[:,sys.argv[3]]
Y= Y
Yerr= data[:,sys.argv[4]]
dX= X[5]-X[4]

print X, Y, Yerr

# Begin of fit code.
# giving initial parameters
sigma = 200
height = 30
p= [sigma,height]

# define your function:
def f(x, s, h ): return h*x*x*exp(-0.5*(x/s)**2)

# fit! (given that data is an array with the data to fit)
# do fit using Levenberg-Marquardt

popt, punc, rchi2, dof = general_fit( f, X, Y, p, Yerr)


for i in range( len(p) ):
  print "%2i %12f +/- %10f"%(i,popt[i],punc[i])


print popt[0] , "+/-", punc[0] 
# Plot

P= np.arange(0,200,1)
# plot data
# errorbar(freq, vr/v0, yerr=uvr, fmt='r+', label="Data")
# plot fit
pl.errorbar( X, Y, yerr=Yerr, fmt='g-', label="Start")
#pl.plot(X, Y, 'g-', label="Start")
pl.plot(P, f(P, popt[0], popt[1]), 'b-', label="Fit")
#pl.yscale("log")
#xlabel("Frequency [Hz]")
#ylabel("Relative amplitude")
#legend(("data","fit"))
#legend()
pl.show()


