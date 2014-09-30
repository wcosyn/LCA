#!/usr/bin/python

# OUTPUTFILE First_File Second_File ..

import sys
import numpy as np

print sys.argv

max= 126
data1 = np.loadtxt(sys.argv[2])
X1 = data1[:max,0]
Y1 = data1[:max,1]
Z1 = data1[:max,2]

for i in range( 3 , len(sys.argv) ):
  data2 = np.loadtxt(sys.argv[i])
  X2 = data2[:max,0]
  Y2 = data2[:max,1]
  Z2 = data2[:max,2]
  if X1[9] != X2[9] :
    print "ERROR X1 != X2"
    exit()
  if Y1[24] != Y2[24] :
    print "ERROR Y1 != Y2"
    exit()
  Z1+= Z2

file = open( sys.argv[1], 'w' )
integral_mf= 0
integral_corr= 0
dX = X1[1]-X1[0]
for i in range(0, len(X1) ):
  file.write( "%f \t\t %f \t\t %f \n" % (X1[i], Y1[i], Z1[i] ) )
  integral_mf+= X1[i]*X1[i]*dX* Y1[i]
  integral_corr+= X1[i]*X1[i]*dX* Z1[i]

file.close()

print "Integral MF : %f" %( integral_mf )
print "Integral CORR : %f" %( integral_corr )
print "Integral CORR : %f" %( integral_mf+ integral_corr )

 
