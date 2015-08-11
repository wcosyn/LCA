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

data_p= np.loadtxt("dens_ob3.00.110.%s.p"%(Sy))
X= data_p[:41,0]
dX= X[1]-X[0]
p_mf= data_p[:41,1]
p_co= data_p[:41,2]
p= p_mf+ p_co
p_check=dX*sum(X*X*(p))
print "norm of p ", p_check

data_n= np.loadtxt("dens_ob3.00.110.%s.n"%(Sy))
X= data_n[:41,0]
dX= X[1]-X[0]
p_mf= data_p[:41,1]
n_mf= data_n[:41,1]
n_co= data_n[:41,2]
n= n_mf+ n_co
n_check=dX*sum(X*X*(n))
print "norm of n ", n_check

print "X range " , X[0], X[-1]

kin_p_mf=dX*sum(X*X*X*X*p_mf)
kin_p_co=dX*sum(X*X*X*X*p_co)
kin_n_mf=dX*sum(X*X*X*X*n_mf)
kin_n_co=dX*sum(X*X*X*X*n_co)
print "P (MeV)\t\t", kin_p_mf*20.7498/Z , "\t", kin_p_co*20.7498/Z, "\t", (kin_p_mf+ kin_p_co)*20.7498/Z
print "N (MeV)\t\t", kin_n_mf*20.7213/N , "\t", kin_n_co*20.7213/N, "\t", (kin_n_mf+ kin_n_co)*20.7213/N
print "Z/N\t\t", kin_p_mf*20.7498/Z/kin_n_mf/20.7213*N, "\t", (kin_p_mf+kin_p_co)*20.7498/Z/(kin_n_mf+kin_n_co)/20.7213*N
print float(Z)/A



print "\n\n\n"
print "N (MeV)\t\t", kin_n_mf*20.7213/N , "\t", kin_n_co*20.7213/N, "\t", (kin_n_mf+ kin_n_co)*20.7213/N
print "P  w int (MeV)\t", kin_p_mf*20.7498/p_check , "\t", kin_p_co*20.7498/p_check, "\t", (kin_p_mf+ kin_p_co)*20.7498/p_check
print "N w int (MeV)\t", kin_n_mf*20.7213/n_check , "\t", kin_n_co*20.7213/n_check, "\t", (kin_n_mf+ kin_n_co)*20.7213/n_check

print "Z/N w int\t", kin_p_mf*20.7498/p_check/kin_n_mf/20.7213*n_check, "\t", (kin_p_mf+kin_p_co)*20.7498/p_check/(kin_n_mf+kin_n_co)/20.7213*n_check


