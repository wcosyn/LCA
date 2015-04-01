#!/usr/bin/python
import sys
import numpy as np
import pylab as pl
import matplotlib.cm as cm
from matplotlib import rc
rc('text', usetex=True)
rc('font', size=24)
#rc('font',family='serif')

len= len(sys.argv)-1 # WUT? seriously Maarten?!?
fig=pl.figure()
fig.subplots_adjust(bottom=0.15, left=0.15)

data = np.loadtxt( sys.argv[2], usecols=[0,1] )
X = data[:,0]
previous = data[:,1]
pl.fill_between( X*np.sqrt(2), 0, previous, color=cm.hot(1./(len) ))
pl.plot(X*np.sqrt(2), previous, color=cm.hot(1./(len) ), label=r"$l=0$" )

for i in range( 3, len ):
  print i
  print sys.argv[i]
  new = np.loadtxt( sys.argv[i], usecols=[1] )
  next = previous+ new
  pl.fill_between( X*np.sqrt(2), previous, next, color=cm.hot(float(i-1)/(len) ), label="l=%d" % (i-2) )
  pl.plot(X*np.sqrt(2), next, color=cm.hot( float(i-1)/(len) ), label=r"$l\leq%d$" % (i-2) )
  previous = next

total = np.loadtxt( sys.argv[len], usecols=[1] )
pl.plot( X*np.sqrt(2), total , color="black", label=r"total")
pl.legend()
#pl.xlabel(r" $P =  \frac{1}{\sqrt{2}} ( p_1+p_2 )$")
pl.xlabel(r" $P =  ( p_1+p_2 )$")
pl.ylabel( "$P_2$" )
normalization = sum( total*X*X) * (X[1]-X[0] )
print "Normalization is ", normalization

pl.show()
#pl.savefig( sys.argv[1] )




