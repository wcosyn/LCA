#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import glob
import os

A="Be"

test_data= np.loadtxt("dens_ob3.00.110.%s.0"%(A))
X= test_data[:,0]
newY= np.zeros_like(X)
newZ= np.zeros_like(X)

newf= open("dens_ob3.00.110.{:s}.-1".format(A), "w")

for files in glob.glob("dens_ob3.00.110.{:s}.?".format(A)):
    data=np.loadtxt(files)
    Y= data[:,1]
    Z= data[:,2]
    newY+= Y
    newZ+= Z

for i in range(0, len(X)):
  newf.write("%f\t%f\t%f\n" % ( X[i], newY[i], newZ[i] ) )

newf.close()
dX= X[1]-X[0]
integral=dX*sum(X*X*(newY+newZ))
print integral


