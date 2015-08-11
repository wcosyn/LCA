import numpy as np
import pylab as pl
import scipy.integrate


nuc = "C"

filename_Y="./{:s}_relative.txt".format(nuc) 
MF_Y =  np.loadtxt(filename_Y,unpack=False)[:,1]
     

filename_X = "../relative/{:s}_tb_rel_new.txt".format(nuc) 
X = np.loadtxt(filename_X,unpack=False)
     
     #fig = pl.figure(figsize=(6,6))
     #fig.subplots_adjust(bottom=0.22,left=0.22)
     #ax  = fig.add_subplot(111)
     
     

k_X =  X[0:200,0]
MF_X = X[0:200,1]
MF_Y = MF_Y[0:200]
dist_Y = MF_Y*k_X*k_X
dist_X = MF_X*k_X*k_X
norm_Y = scipy.integrate.trapz(dist_Y,x=None,dx = 0.02)
norm_X = scipy.integrate.trapz(dist_X,x=None,dx = 0.02)
print (MF_Y)/(MF_X/norm_X)
print norm_Y
print norm_X
