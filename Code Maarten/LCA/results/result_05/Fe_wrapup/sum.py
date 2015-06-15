#!/usr/bin/python

import sys
import numpy as np
import pylab as pl

data= np.loadtxt("dens_rel_all.E2.Fe3", usecols=[0,1] )
X= data[:,0]
Y= data[:,1]
print len(Y), Y[0]

for i in range( 3, 7 ):
  print i
  data= np.loadtxt("dens_rel_all.E%d.Fe3" % (i), usecols=[1] )
  print len(data), data[0]+ Y[0]
  Y+= data
  print Y[0]

f = open( "dens_rel_all.E2+.Fe3", "w" )
f.write( "momentum [fm] \t relative momentum distribution n_2^(E), E>=2 \n")
for i in range(0, len(X) ):
  f.write("%f\t%f\n" % ( X[i], Y[i] ) )
  print X[i], Y[i]

f.close()

