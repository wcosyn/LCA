#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import glob
import os

A="C"

test_data= np.loadtxt("dens_ob3.00.110.%s.0"%(A))
X= test_data[:,0]
newY= np.zeros_like(X)
newZ= np.zeros_like(X)

#newf= open("dens_ob3.00.110.{:s}.-1".format(A), "w")

for files in glob.glob("dens_ob3.00.110.{:s}.?".format(A)):
    data=np.loadtxt(files)
    Y= data[:,1]
    Z= data[:,2]
    newY+= Y
    newZ+= Z

data=np.loadtxt("dens_ob3.00.110.{:s}.-3".format(A))
x_all= data[16:40,0]
all= data[16:40,1]+ data[16:40,2]

data=np.loadtxt("c12.momentum.txt")
x_argonne= data[16:40,0]
argonne= data[16:40,1]*2/19.7392088


dX= X[1]-X[0]
integral=dX*sum(X*X*(newY+newZ))
print "0-9 sum integral", integral
dX= x_all[1]-x_all[0]
integral_all=dX*sum(x_all*x_all*(all))
print "all integral", integral_all
dX= x_argonne[1]-x_argonne[0]
integral_argonne=  dX* sum( x_argonne* x_argonne* argonne )
print "argonne int", integral_argonne



