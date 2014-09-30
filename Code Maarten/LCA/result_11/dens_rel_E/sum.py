#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import os.path

data= np.loadtxt("dens_rel.00.111.%s.3" %(sys.argv[1]), usecols=[0,3] )
X= data[:,0]
Y= data[:,1]
print len(Y), Y[0]

for i in range( 4, 8 ):
  file="dens_rel.00.111.{:s}.{:d}".format(sys.argv[1], i)
  if( os.path.exists(file)):
    data= np.loadtxt(file, usecols=[3] )
    print len(data), data[0]+ Y[0]
    Y+= data
    print Y[0]

f = open( "dens_rel.00.111.%s.3+" %(sys.argv[1]) , "w" )
f.write( "momentum [fm] \t relative momentum distribution n_2^(E), E>=3 \n")
for i in range(0, len(X) ):
  f.write("%f\t%f\n" % ( X[i], Y[i] ) )
#  print X[i], Y[i]

f.close()

