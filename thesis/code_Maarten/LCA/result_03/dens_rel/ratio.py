#!/usr/bin/python

import io
import sys
import math as m
import numpy as np
import pylab as pl

Y1=[]
Y2=[]
X=[]

data1= np.loadtxt( sys.argv[1], usecols=[0,1,2] )
data2= np.loadtxt( sys.argv[2], usecols=[0,1,2] )

Y1= data1[:,1]
Y2= data2[:,1]
Z1= data1[:,2]
Z2= data2[:,2]
X= data1[:,0]

pl.plot( X, Y1/Y2 )
pl.plot( X, Z1/Z2 )


pl.show()
