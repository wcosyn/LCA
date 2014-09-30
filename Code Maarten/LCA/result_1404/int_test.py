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


for files in glob.glob("{:s}/dens_ob3.00.110.{:s}.?".format(A,A)):
    print files
    data=np.loadtxt(files)
    mf= data[:,1]
    co= data[:,2]
    Z= data[:,3]
    summf+= mf
    sumco+= co
    sumZ+= Z

dX= X[1]-X[0]
integral=dX*sum(X*X*(sumZ))
integralmf=dX*sum(X*X*(summf))
integralco=dX*sum(X*X*(sumco))

print "TEST sum integral: ", integral
print "TEST sum mf integral: ", integralmf
print "TEST sum co integral: ", integralco

data_all= np.loadtxt( "{:s}/dens_ob3.00.110.{:s}.-1".format(A,A) )
mf = data_all[:,1]
co = data_all[:,2]
Z = data_all[:,3]
to_integralmf=dX*sum(X*X*(mf))
to_integralco=dX*sum(X*X*(co))
to_integral=dX*sum(X*X*(Z))
print "TEST total integral: ", to_integral, "/", integral, ":\t", to_integral/integral
print "TEST total mf integral: ", to_integralmf, "/", integralmf, ":", to_integralmf/integralmf
print "TEST total co integral: ", to_integralco, "/", integralco, ":", to_integralco/integralco


