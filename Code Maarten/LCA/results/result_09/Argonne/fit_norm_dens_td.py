#!/usr/bin/python
import sys
import pylab as pl
import numpy as np
from scipy.optimize import curve_fit




# Argonne data
column= int(sys.argv[4])
xa, ya, dya= np.loadtxt("table."+sys.argv[2], usecols=(0,column,column+1), unpack=True)
ya=ya/389.6363
dya=dya/389.6363

# My data
P_all, k_all, y_all= np.loadtxt(sys.argv[3], usecols=(0,1,4), unpack=True)

k= []
y= []
P= float(sys.argv[1])
for i in range(0, len(P_all)):
  if (P_all[i] == P ):
    k.append(k_all[i])
    y.append(y_all[i])
dk = k[1]-k[0]

def get_k( x ):
  intx= int(np.floor(x/dk))
  y0= y[intx]
  if( intx < len(y)-1 ):
    y1= y[intx+1]
  else:
    y1=0
  dy= y1-y0
  dx= x- intx*dk
  return (y0+ dy*dx/dk)

k_select=[]
y_select=[]
for i in range(0, len(xa)):
  k_select.append(xa[i])
  y_select.append(get_k(xa[i]))

#print k[0:20]
#print y[0:20]
#print k_select[0:10]
#print y_select[0:10]

  
def func( x,a ):
  return [a*y for y in y_select]

popt, pcov= curve_fit( func, xa, ya, p0=1.1)
sys.stdout.write(`popt[0]`)

#print popt
#print pcov
#pl.plot( xa, ya, 'ro', xa,func(k,popt), 'b-')
#pl.yscale('log')
#pl.show()








