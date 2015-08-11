#!/usr/bin/python
# Different with simulate.py is that ws calculations give P=p1+p2
# while ho calculations give P = (p1+p2)/sqrt(2)

import io
import sys
import math as m
import numpy as np
import pylab as pl


def P2( x ):
  i = int( x/dx )
  if( i >= len(F) ):
    return 0
  return F[i]

def draw_p2():
  while True:
    r = 1*np.random.random_sample()
    f= P2( r )
    distr = max*1
    frac= f/ distr
    if frac > 1:
      print "ERROR", r, P2(r), distr
      continue
    y = np.random.random_sample()
    if frac  >  y:
      return r
      
def draw_cosTh( ):
  return 1-2*np.random.random_sample()

def draw_phi():
  return 2*m.pi*np.random.random_sample()

x= []
y= []
z= []
p_list= []
data= np.loadtxt( sys.argv[1], usecols= [0,int(sys.argv[2])] )
print sys.argv[1]
sys.stdout.flush()
X= data[:,0]
#X=X*m.sqrt(2)
Y= data[:,1]
#Y=Y/m.sqrt(8)
F=X*X*Y
dx= X[1]-X[0]
norm= sum(F)*dx
F /= norm
max= np.amax( F )
print "norm ",  norm


for i in range(0,100000):
  p =  draw_p2()
  p_list.append( p )
  u= draw_cosTh()
  phi= draw_phi()

  x.append( p* m.sqrt(1-u*u)* m.cos( phi ) )
  y.append( p* m.sqrt(1-u*u)* m.sin( phi ) )
  z.append( p* u )

n, bins, patches= pl.hist( p_list, 20, normed=True )
pl.plot( X, F  )
#pl.show()

mean= sum(X*F)/ sum(F) # Mean
moment2= sum( np.power(X-mean, 2)*F) / sum(F) # Variance
std_dev = np.sqrt(moment2)
skewness= sum( np.power( (X-mean), 3)* F)/ sum(F)/ m.pow( std_dev, 3) # skewness
kurtosis= sum( np.power( (X-mean), 4)* F)/ sum(F)/ m.pow( std_dev, 4) - 3 # kurtosis

print "P2 mean", mean
print "P2 std_dev", std_dev
print "P2 skewness", skewness
print "P2 kurtosis", kurtosis



F, bins, patches = pl.hist( x, 20, normed=True )
#pl.show()
x_values= []
for i in range(1, len(bins) ):
  x_values.append( 0.5*( bins[i-1]+bins[i] ) )

mean= sum(x_values*F)/ sum(F) # Mean
moment2= sum( np.power(x_values-mean, 2)*F) / sum(F) # Variance
std_dev = np.sqrt(moment2)
skewness= sum( np.power( (x_values-mean), 3)* F)/ sum(F)/ m.pow( std_dev, 3) # skewness
kurtosis= sum( np.power( (x_values-mean), 4)* F)/ sum(F)/ m.pow( std_dev, 4) - 3 # kurtosis
print "x mean", mean
print "x std_dev", std_dev
print "x skewness", skewness
print "x kurtosis", kurtosis





F, bins, patches = pl.hist( y, 20, normed=True )
#pl.show()
x_values= []
for i in range(1, len(bins) ):
  x_values.append( 0.5*( bins[i-1]+bins[i] ) )

mean= sum(x_values*F)/ sum(F) # Mean
moment2= sum( np.power(x_values-mean, 2)*F) / sum(F) # Variance
std_dev = np.sqrt(moment2)
skewness= sum( np.power( (x_values-mean), 3)* F)/ sum(F)/ m.pow( std_dev, 3) # skewness
kurtosis= sum( np.power( (x_values-mean), 4)* F)/ sum(F)/ m.pow( std_dev, 4) - 3 # kurtosis
print "y mean", mean
print "y std_dev", std_dev
print "y skewness", skewness
print "y kurtosis", kurtosis



F, bins, patches = pl.hist( z, 20, normed=True )
#pl.show()
x_values= []
for i in range(1, len(bins) ):
  x_values.append( 0.5*( bins[i-1]+bins[i] ) )

mean= sum(x_values*F)/ sum(F) # Mean
moment2= sum( np.power(x_values-mean, 2)*F) / sum(F) # Variance
std_dev = np.sqrt(moment2)
skewness= sum( np.power( (x_values-mean), 3)* F)/ sum(F)/ m.pow( std_dev, 3) # skewness
kurtosis= sum( np.power( (x_values-mean), 4)* F)/ sum(F)/ m.pow( std_dev, 4) - 3 # kurtosis
print "z mean", mean
print "z std_dev", std_dev
print "z skewness", skewness
print "z kurtosis", kurtosis





#print "means ", meanx, meany, meanz
#print "std_dev", std_devx, std_devy, std_devz
#print "skewness, ", skewnessx, skewnessy, skewnessz
#print "kurtosis, ", kurtosisx, kurtosisy, kurtosisz

#x= np.arange(0, 3, 0.01)
#pl.plot( x, MB(x,std_dev) )
    





