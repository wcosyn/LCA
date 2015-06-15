#!/usr/bin/python
#
# Usage: ./kinenergy 9 4 Be
#

import sys
import numpy as np
import pylab as pl
import glob
import os


A= int(sys.argv[1])
Z= int(sys.argv[2])
Sy= sys.argv[3]
N= A-Z

print Sy

for line in open("{:s}/dens_ob4.00.111.{:s}.-1-1-1-1.p".format(Sy, Sy)):
  if "norm" in line:
    norm= float(line.split("=")[-1])
    break

data_p= np.loadtxt("{:s}/dens_ob4.00.111.{:s}.-1-1-1-1.p".format(Sy, Sy))
X= data_p[:41,0]
dX= X[1]-X[0]
p_mf= data_p[:41,1]
p_co= data_p[:41,2]

p= p_mf+ p_co
p_check=dX*sum(X*X*(p))
p_check_mf=dX*sum(X*X*(p_mf))
print "norm of p ", p_check
print "mf norm of p ", p_check_mf*norm

data_n= np.loadtxt("{:s}/dens_ob4.00.111.{:s}.-1-1-1-1.n".format(Sy, Sy))
X= data_n[:41,0]
dX= X[1]-X[0]
n_mf= data_n[:41,1]
n_co= data_n[:41,2]

n= n_mf+ n_co
n_check=dX*sum(X*X*(n))
n_check_mf=dX*sum(X*X*(n_mf))
print "norm of n ", n_check
print "mf norm of n ", n_check_mf*norm

print "X range " , X[0], X[-1]

p_mf= data_p[:31,1]
n_mf= data_n[:31,1]
X_mf= data_p[:31,0]

kin_p_mf=dX*sum(X_mf*X_mf*X_mf*X_mf*p_mf)
kin_p_co=dX*sum(X*X*X*X*p_co)
kin_n_mf=dX*sum(X_mf*X_mf*X_mf*X_mf*n_mf)
kin_n_co=dX*sum(X*X*X*X*n_co)
print "P (MeV)\t\t", kin_p_mf*20.7498/Z*norm , "\t", kin_p_co*20.7498/Z, "\t", (kin_p_mf+ kin_p_co)*20.7498/Z
print "N (MeV)\t\t", kin_n_mf*20.7213/N*norm , "\t", kin_n_co*20.7213/N, "\t", (kin_n_mf+ kin_n_co)*20.7213/N
print "Z/N\t\t", kin_p_mf*20.7498/Z/kin_n_mf/20.7213*N, "\t", (kin_p_mf+kin_p_co)*20.7498/Z/(kin_n_mf+kin_n_co)/20.7213*N
print float(Z)/A




