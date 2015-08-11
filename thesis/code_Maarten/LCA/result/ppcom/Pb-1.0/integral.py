#!/usr/bin/python

import sys
import numpy as np


data1 = np.loadtxt(sys.argv[1])
#data2 = np.loadtxt(sys.argv[2])

X1 = data1[:,0]
Y1 = data1[:,1]
Y2 = data1[:,2]

#if X1[9] != X2[9] :
#	print "ERROR X1 != X2"
#	exit()

dX= X1[2]-X1[1]

integral1= 0
integral2= 0
for i in range(0, len(X1) ):
	integral1+= X1[i]*X1[i]*dX* Y1[i]
	integral2+= X1[i]*X1[i]*dX* Y2[i]

print " %f\t%s" %( integral1, sys.argv[1])

