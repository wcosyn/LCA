#!/usr/bin/python

# New_to_add_file old_file sum_file factor_for_new_to_add_file

import sys
import numpy as np



data1 = np.loadtxt(sys.argv[1])
data2 = np.loadtxt(sys.argv[2])

file = open( sys.argv[3], 'w' )

X1 = data1[:,0]
Y1 = data1[:,1]
Y12 = data1[:,2]
X2 = data2[:,0]
Y2 = data2[:,1]
Y22 = data2[:,2]

factor= float( sys.argv[4] )

if X1[9] != X2[9] :
	print "ERROR X1 != X2"
	exit()

dX= X1[2]-X2[1]
newY = factor*Y1+Y2
newY2 = factor*Y12+Y22
integral= 0
for i in range(0, len(X1) ):
	file.write( "%f \t\t %f \t\t %f \n" % (X1[i], newY[i], newY2[i] ) )
	integral+= X2[i]*X2[i]*dX* newY[i]

file.close()

print "Integral: %f" %( integral )

 
