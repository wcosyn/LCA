#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import glob
import os

A="Ag"

test_data= np.loadtxt("dens_ob3.11.110.{0}.0".format(A), usecols=(1,))
mf= np.zeros_like( test_data )
co= np.zeros_like( test_data )
X= np.zeros_like( test_data )
mf3p= np.zeros_like( test_data )
co3p= np.zeros_like( test_data )

int_all= 0
for i in np.arange( 0, 9, 1):
  X= np.loadtxt("dens_ob3.11.110.{0}.{1}".format(A,i), usecols=(0,))
  pp_data= np.loadtxt("dens_ob3.11.110.{0}.{1}".format(A,i), usecols=(1,2))
  nn_data= np.loadtxt("dens_ob3.-1-1.110.{0}.{1}".format(A,i), usecols=(1,2))
  pn_data= np.loadtxt("dens_ob3.-11.110.{0}.{1}".format(A,i), usecols=(1,2))
  all_data= pp_data+nn_data+pn_data
  mf+= all_data[:,0]
  co+= all_data[:,1]
  if i >= 3:
    mf3p+= all_data[:,0]
    co3p+= all_data[:,1]

  newf= open("dens_ob3.00.110.{0}.{1}.b".format(A,i), "w")
  for j in np.arange( 0, len(all_data[:,0])):
    newf.write("{0}\t{1}\t{2}\n".format(X[j], all_data[j,0], all_data[j,1]) )

  newf.close()
  dX= X[1]-X[0]
  integral=dX*sum(X*X*(all_data[:,0]+ all_data[:,1]))
  print i, integral
  int_all+= integral


newf_all= open("dens_ob3.00.110.{0}.-3.b".format(A), "w")
for j in np.arange( 0, len(X)):
  newf_all.write("{0}\t{1}\t{2}\n".format(X[j], mf[j], co[j]) )

newf3p= open("dens_ob3.00.110.{0}.3+.b".format(A), "w")
for j in np.arange( 0, len(X)):
  newf3p.write("{0}\t{1}\t{2}\n".format(X[j], mf3p[j], co3p[j]) )

newf3p.close()
newf_all.close()


print "total int", int_all




