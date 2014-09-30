#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import glob
import os

A=sys.argv[1]

data_all= np.loadtxt("./{:s}/dens_ob4.00.111.{:s}.-1-1-1-1".format(A,A))
X= data_all[:,0]
all_mf= data_all[:,1]
all_co= data_all[:,2]
all= data_all[:,3]

data_pp= np.loadtxt("./{:s}/dens_ob4.11.111.{:s}.-1-1-1-1".format(A,A))
data_nn= np.loadtxt("./{:s}/dens_ob4.-1-1.111.{:s}.-1-1-1-1".format(A,A))
data_pn= np.loadtxt("./{:s}/dens_ob4.-11.111.{:s}.-1-1-1-1".format(A,A))
pp_mf= data_pp[:,1]
pp_co= data_pp[:,2]
pp_all= data_pp[:,3]
nn_mf= data_nn[:,1]
nn_co= data_nn[:,2]
nn_all= data_nn[:,3]
pn_mf= data_pn[:,1]
pn_co= data_pn[:,2]
pn_all= data_pn[:,3]

pp_ra= (pp_all)/all
nn_ra= (nn_all)/all
pn_ra= (pn_all)/all

newf= open("./dens_ratio.111.%s.-1-1-1-1"%(A), "w")
for i in range(0, len(X)):
  newf.write("%f\t%f\t%f\t%f\n" % ( X[i], pp_ra[i], nn_ra[i], pn_ra[i] ) )

newf.close()
dX= X[1]-X[0]
integral=dX*sum(X*X*(pp_ra+nn_ra+pn_ra))
integral=dX*sum(X*X*(pp_ra+nn_ra+pn_ra)*all)
print integral


