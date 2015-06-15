#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import glob
import os

A=sys.argv[1]

data_all= np.loadtxt("dens_ob3.00.110.%s.-1"%(A))
X= data_all[:,0]
all_mf= data_all[:,1]
all_co= data_all[:,2]
all= all_mf+all_co

data_pp= np.loadtxt("dens_ob3.11.110.%s.-1"%(A))
data_nn= np.loadtxt("dens_ob3.-1-1.110.%s.-1"%(A))
data_pn= np.loadtxt("dens_ob3.-11.110.%s.-1"%(A))
pp_mf= data_pp[:,1]
pp_co= data_pp[:,2]
nn_mf= data_nn[:,1]
nn_co= data_nn[:,2]
pn_mf= data_pn[:,1]
pn_co= data_pn[:,2]

pp_ra= (pp_mf+pp_co)/all
nn_ra= (nn_mf+nn_co)/all
pn_ra= (pn_mf+pn_co)/all

newf= open("dens_ratio.110.%s.-1"%(A), "w")
for i in range(0, len(X)):
  newf.write("%f\t%f\t%f\t%f\n" % ( X[i], pp_ra[i], nn_ra[i], pn_ra[i] ) )

newf.close()
dX= X[1]-X[0]
integral=dX*sum(X*X*(pp_ra+nn_ra+pn_ra))
integral=dX*sum(X*X*(pp_ra+nn_ra+pn_ra)*all)
print integral


