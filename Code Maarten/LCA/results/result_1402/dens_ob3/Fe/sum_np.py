#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import glob
import os

A="Fe"



X= np.loadtxt("dens_ob3.11.110.{0}.{1}".format(A,-1), usecols=(0,))
pp_data= np.loadtxt("dens_ob3.11.110.{0}.{1}".format(A,-1), usecols=(1,2))
nn_data= np.loadtxt("dens_ob3.-1-1.110.{0}.{1}".format(A,-1), usecols=(1,2))
pn_data= np.loadtxt("dens_ob3.-11.110.{0}.{1}".format(A,-1), usecols=(1,2))

mf_n= nn_data[:,0]+ 0.5*pn_data[:,0]
co_n= nn_data[:,1]+ 0.5*pn_data[:,1]
mf_p= pp_data[:,0]+ 0.5*pn_data[:,0]
co_p= pp_data[:,1]+ 0.5*pn_data[:,1]

dX= X[1]-X[0]


newf_n= open("dens_ob3.00.110.{0}.n.b".format(A), "w")
for j in np.arange( 0, len(X)):
  newf_n.write("{0}\t{1}\t{2}\n".format(X[j], mf_n[j], co_n[j]) )

newf_p= open("dens_ob3.00.110.{0}.p.b".format(A), "w")
for j in np.arange( 0, len(X)):
  newf_p.write("{0}\t{1}\t{2}\n".format(X[j], mf_p[j], co_p[j]) )

newf_p.close()
newf_n.close()

int_n= dX*sum(X*X*(mf_n+ co_n))
int_p= dX*sum(X*X*(mf_p+ co_p))

print "norm n", int_n
print "norm p", int_p
print "norm total", int_n+ int_p




