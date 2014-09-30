#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import glob
import os

A=sys.argv[1]

test_data= np.loadtxt("%s/dens_ob3.00.110.%s.0"%(A,A))
X= test_data[:,0]
summf= np.zeros_like(X)
sumco= np.zeros_like(X)
sumZ= np.zeros_like(X)

#newf= open("{:s}/dens_ob3.00.110.{:s}.-1".format(A,A), "w")
newf= open("./dens_ob3.00.110.{:s}.-1".format(A), "w")

for files in glob.glob("{:s}/dens_ob3.00.110.{:s}.?".format(A,A)):
    print files
    data=np.loadtxt(files)
    mf= data[:,1]
    co= data[:,2]
    Z= data[:,3]
    summf+= mf
    sumco+= co
    sumZ+= Z


for i in range(0, len(X)):
  newf.write("%f\t%f\t%f\t%f\n" % ( X[i], summf[i], sumco[i], sumZ[i] ) )

dX= X[1]-X[0]
integral=dX*sum(X*X*(sumZ))
print "TEST sum integral: ", integral

#
# IDEM FOR n
#

test_data= np.loadtxt("%s/dens_ob3.00.110.%s.0.n"%(A,A))
X= test_data[:,0]
summf= np.zeros_like(X)
sumco= np.zeros_like(X)
sumZ= np.zeros_like(X)

#newf= open("{:s}/dens_ob3.00.110.{:s}.-1".format(A,A), "w")
newf= open("./dens_ob3.00.110.{:s}.-1.n".format(A), "w")

for files in glob.glob("{:s}/dens_ob3.00.110.{:s}.?.n".format(A,A)):
    print files
    data=np.loadtxt(files)
    mf= data[:,1]
    co= data[:,2]
    Z= data[:,3]
    summf+= mf
    sumco+= co
    sumZ+= Z


for i in range(0, len(X)):
  newf.write("%f\t%f\t%f\t%f\n" % ( X[i], summf[i], sumco[i], sumZ[i] ) )

dX= X[1]-X[0]
nintegral=dX*sum(X*X*(sumZ))
print "TEST sum n integral: ", nintegral

#
# IDEM FOR p
#

test_data= np.loadtxt("%s/dens_ob3.00.110.%s.0.p"%(A,A))
X= test_data[:,0]
summf= np.zeros_like(X)
sumco= np.zeros_like(X)
sumZ= np.zeros_like(X)

#newf= open("{:s}/dens_ob3.00.110.{:s}.-1".format(A,A), "w")
newf= open("./dens_ob3.00.110.{:s}.-1.p".format(A), "w")

for files in glob.glob("{:s}/dens_ob3.00.110.{:s}.?.p".format(A,A)):
    print files
    data=np.loadtxt(files)
    mf= data[:,1]
    co= data[:,2]
    Z= data[:,3]
    summf+= mf
    sumco+= co
    sumZ+= Z


for i in range(0, len(X)):
  newf.write("%f\t%f\t%f\t%f\n" % ( X[i], summf[i], sumco[i], sumZ[i] ) )

dX= X[1]-X[0]
pintegral=dX*sum(X*X*(sumZ))
print "TEST sum p integral: ", pintegral

print nintegral+pintegral, integral


