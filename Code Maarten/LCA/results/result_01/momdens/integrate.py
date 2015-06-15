#!/usr/bin/python
import sys
import numpy as np
import pylab as pl

data= np.loadtxt(sys.argv[1])
r= data[:,0]
#r= r*np.sqrt(2)
dr = np.zeros( [len(r) ] )
#dr[0] = r[0]
#for i in range(1, len(r) ):
#  dr[i] = r[i] - r[i-1]

for i in range( 0, len(r)-1 ):
  dr[i] = r[i+1] - r[i]


print dr


for i in range(1,len(data[0])):
  y = data[:,i]
#  y = y/sqrt(8)
#  integral= sum( y*r*r )*dr
  integral= sum( y *dr)
  print "The integral over column", i, "is", integral


