#!/usr/bin/python
# ARGS: INPUTFILE
import sys
import io
import math as m
import numpy as np


# Read column from file
data = np.loadtxt(sys.argv[1], usecols= [0,1]  )

X = data[:,0] 
# From computation units to experimental units
X = X*np.sqrt(2)

Y = data[:,1]
F = X*X*Y

#Calculate moments

mean= sum(X*F)/ sum(F) # Mean
moment2= sum( np.power(X-mean, 2)*F) / sum(F) # Variance
std_dev = np.sqrt(moment2)
#moment2= sum( np.power(X, 2)*F) / sum(F) # raw 2nd moment 
#std_dev= m.sqrt(moment2 - moment1*moment1 )

skewness= sum( np.power( (X-mean), 3)* F)/ sum(F)/ m.pow( std_dev, 3) # skewness
kurtosis= sum( np.power( (X-mean), 4)* F)/ sum(F)/ m.pow( std_dev, 4) - 3 # kurtosis

#a1= moment1/2.*m.sqrt(m.pi/2.)
#a2= std_dev*m.sqrt(m.pi/(3*m.pi-8) )
#print "a1, a2 ", a1, a2
print "mean", mean
print "2nd moment", moment2
print "std_dev", std_dev
print "Px std_dev", std_dev*m.sqrt(m.pi/(3*m.pi-8))
print "skewness", skewness
print "kurtosis", kurtosis

lam= mean*mean*mean/moment2
alpha= mean*mean/moment2
beta= mean/moment2
print "mu, lambda ", mean, lam
print "alpha, beta ", alpha, beta 
print
print sum(F ) * (X[1]-X[0])/m.sqrt(8)

 

