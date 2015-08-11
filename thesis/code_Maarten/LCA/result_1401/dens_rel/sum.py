#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import os.path

data= np.loadtxt("sorted_dens_rel.00.110.%s.0" %(sys.argv[1]), usecols=[0,4] )
X= data[:,0]
Y= data[:,1]

Ytp= np.zeros_like( Y )

for i in range( 1, 8 ):
  file="sorted_dens_rel.00.110.{:s}.{:d}".format(sys.argv[1], i)
  if( os.path.exists(file)):
    data= np.loadtxt(file, usecols=[4] )
    print len(data), data[0]+ Y[0]
    Y+= data
    if i >= 3:
      Ytp += data

if Ytp[3] != 0:
  f = open( "sorted_dens_rel.00.110.%s.3+" %(sys.argv[1]) , "w" )
  f.write( "momentum [fm] \t relative momentum distribution n_2^(E), E>=3 \n")
  for i in range(0, len(X) ):
    f.write("%f\t%f\n" % ( X[i], Ytp[i] ) )


f = open( "sorted_dens_rel.00.110.%s.-1" %(sys.argv[1]) , "w" )
f.write( "momentum [fm] \t relative momentum distribution n_2^(All E) \n")
for i in range(0, len(X) ):
  f.write("%f\t%f\n" % ( X[i], Y[i] ) )

f.close()

