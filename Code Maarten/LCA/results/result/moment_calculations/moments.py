#!/usr/bin/python
# ARGS: INPUTFILE
import sys
import io
import math as m
import numpy as np


# Read column from file
data = np.loadtxt(sys.argv[1], usecols= [0,1]  )

Xp = data[:,0] 
#print Xp 
Xm = -1.*Xp
Xr = Xm[::-1].copy()
s = Xr.size
X = np.zeros( 2*s-1 )
X[:s-1] = Xr[:s-1]
X[s-1:] = Xp
#print X

Yp = data[:,1]
#print Yp
Ym = Yp
Yr = Ym[::-1].copy()
s = Yr.size
Y = np.zeros( 2*s-1 )
Y[:s-1] = Yr[:s-1]
Y[s-1:] = Yp
F = np.power(Y,1./3.)
#print Y

X = X*np.sqrt(2)

#Calculate moments

moment1= sum( X*Y)/ sum(Y) # Mean
moment2= sum( np.power( (X-moment1), 2)* Y)/ sum(Y) # Variance
std_dev= m.sqrt( moment2 )
#moment2= sum( np.power( (X), 2)* Y)/ sum(Y) # Variance
#std_dev= m.sqrt(moment2 - moment1*moment1 )

skewness= sum( np.power( (X-moment1), 3)* Y)/ sum(Y)/ m.pow( std_dev, 3) # skewness
kurtosis= sum( np.power( (X-moment1), 4)* Y)/ sum(Y)/ m.pow( std_dev, 4) - 3 # kurtosis

print "mean", moment1
print "2nd moment", moment2
print "std_dev", std_dev
print "P2 std. dev", std_dev* m.sqrt( (3.*m.pi-8)/m.pi) 
print "skewness", skewness
print "kurtosis", kurtosis

 

