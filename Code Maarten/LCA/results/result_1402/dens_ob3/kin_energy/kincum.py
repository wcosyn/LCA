#!/usr/bin/python
#
# Usage: ./kinenergy 9 4 Be
#

import sys
import numpy as np
import pylab as pl
import glob
import os


A= int(sys.argv[1])
Z= int(sys.argv[2])
Sy= sys.argv[3]
N= A-Z

print Sy

data_p= np.loadtxt("dens_ob3.00.110.%s.p"%(Sy))

data_n= np.loadtxt("dens_ob3.00.110.%s.n"%(Sy))

newf= open("kin_cum.%s"%(A), "w")
ncum= 0
nmf= 0
ncorr=0
pcum= 0
pmf= 0
pcorr=0

dX = data_n[1,0] - data_n[0,0]

X= data_n[:41,0]
n= data_n[:41,1]+ data_n[:41,2]
p= data_p[:41,1]+ data_p[:41,2]
#norm_p= dX* sum( X*X*X*X*p)*20.7498/Z
#norm_n= dX* sum( X*X*X*X*n)*20.7213/N
norm_p= 1
norm_n= 1

for end in range( 1, 41 ):
  newX= data_n[end,0]
  p_mf= data_p[end,1]
  p_co= data_p[end,2]
  p= p_mf+ p_co
  n_mf= data_n[end,1]
  n_co= data_n[end,2]
  n= n_mf+ n_co

  pcum+= dX*newX*newX*newX*newX*p/norm_p*20.7498/Z
  pmf+= dX*newX*newX*newX*newX*p_mf/norm_p*20.7498/Z
  pcorr+= dX*newX*newX*newX*newX*p_co/norm_p*20.7498/Z
  ncum+= dX*newX*newX*newX*newX*n/norm_n*20.7213/N
  nmf+= dX*newX*newX*newX*newX*n_mf/norm_n*20.7213/N
  ncorr+= dX*newX*newX*newX*newX*n_co/norm_n*20.7213/N

  newf.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format( newX, pcum, pmf,pcorr, ncum, nmf, ncorr ))
