#!/usr/bin/python
import sys
import pylab as pl
import numpy as np
from scipy.optimize import curve_fit




# Argonne data
column= int(sys.argv[4])
xa, ya, dya= np.loadtxt("table."+sys.argv[2], usecols=(0,column,column+1), unpack=True)
dxa=xa[1]-xa[0]


int1= sum(ya*xa*xa)*dxa/2/np.pi/np.pi

# My data
P_all, k_all, y_all= np.loadtxt(sys.argv[3], usecols=(0,1,4), unpack=True)

k= []
y= []
P= float(sys.argv[1])
for i in range(0, len(P_all)):
  if (P_all[i] == P ):
    k.append(k_all[i])
    y.append(y_all[i])

kar = np.array(k)
yar = np.array(y)
dk= kar[1]- kar[0]

int2=sum(yar*kar*kar)*dk*2*np.pi*np.pi


sys.stdout.write(`int1/int2`)

#pl.plot( xa, ya, 'ro', xa,func(k,popt), 'b-')
#pl.yscale('log')
#pl.show()








