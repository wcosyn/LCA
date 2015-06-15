#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import glob
import os

A=sys.argv[1]
print A

test_data= np.loadtxt("./{:s}/dens_ob4.00.110.{:s}.0".format(A,A))
X= test_data[:,0]
new_mf= np.zeros_like(X)
new_co= np.zeros_like(X)
new_all= np.zeros_like(X)

#newf= open("dens_ob3.00.110.{:s}.-1".format(A), "w")

newf_3= open("{:s}/dens_ob4.00.110.{:s}.3+".format(A,A), "w")
new_mf3= np.zeros_like(X)
new_co3= np.zeros_like(X)
new_all3= np.zeros_like(X)

#for files in glob.glob("dens_ob3.00.110.{:s}.?".format(A)):
#    data=np.loadtxt(files)
#    Y= data[:,1]
#    Z= data[:,2]
#    newY+= Y
#    newZ+= Z

for i in range(0, 10):
  fname= "{:s}/dens_ob4.00.110.{:s}.{:d}".format(A,A,i)
  if os.path.isfile(fname):
    print fname
    data=np.loadtxt(fname)
    mf= data[:,1]
    co= data[:,2]
    all= data[:,3]
    new_mf+= mf
    new_co+= co
    new_all+= all
    if i >= 3:
      new_mf3+= mf
      new_co3+= co
      new_all3+= all


for i in range(0, len(X)):
#  newf.write("%f\t%f\t%f\n" % ( X[i], newY[i], newZ[i] ) )
  newf_3.write("%f\t%f\t%f\t%f\n" % ( X[i], new_mf3[i], new_co3[i], new_all3[i] ) )

#newf.close()
newf_3.close()
dX= X[1]-X[0]
integral=dX*sum(X*X*(new_all))
print integral




