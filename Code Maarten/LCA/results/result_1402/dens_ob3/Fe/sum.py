#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import glob
import os

A="Fe"

test_data= np.loadtxt("dens_ob3.00.110.%s.0"%(A))
X= test_data[:,0]
newY= np.zeros_like(X)
newZ= np.zeros_like(X)

newf= open("dens_ob3.00.110.{:s}.-3".format(A), "w")

newf_3= open("dens_ob3.00.110.{:s}.3+".format(A), "w")
newY_3= np.zeros_like(X)
newZ_3= np.zeros_like(X)

#for files in glob.glob("dens_ob3.00.110.{:s}.?".format(A)):
#    data=np.loadtxt(files)
#    Y= data[:,1]
#    Z= data[:,2]
#    newY+= Y
#    newZ+= Z

for i in range(0, 10):
  fname= "dens_ob3.00.110.{:s}.{:d}".format(A,i)
  print fname
  if os.path.isfile(fname):
    data=np.loadtxt(fname)
    Y= data[:,1]
    Z= data[:,2]
    newY+= Y
    newZ+= Z
    if i >= 3:
      newY_3+= Y
      newZ_3+= Z


for i in range(0, len(X)):
  newf.write("%f\t%f\t%f\n" % ( X[i], newY[i], newZ[i] ) )
  newf_3.write("%f\t%f\t%f\n" % ( X[i], newY_3[i], newZ_3[i] ) )

newf.close()
newf_3.close()
dX= X[1]-X[0]
integral=dX*sum(X*X*(newY+newZ))
print integral




