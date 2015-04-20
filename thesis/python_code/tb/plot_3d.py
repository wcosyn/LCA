import numpy as np
import pylab as pl
from scipy import interpolate
import scipy.integrate
import matplotlib.pyplot as plt
import matplotlib.cm as cm
###############################
### global variables   ########
###############################

nuc = ["He","Be","C","O","Al","Ar","Ca40","Fe","Ag","Pb"]
A   = [4,9,12,16,27,40,40,56,109,208]
NORMALISATION_CONVENTION = "A" # choose "A" or "1", yes those are strings!


#################################
# set a bunch of plot parameters#
#################################

params = {'legend.fontsize' : 25,
          'legend.linewidth': 2,
          'axes.linewidth'  : 3.5,
          'axes.labelsize'  : 30,
	  'xtick.minor.size'  : 6,
	  'xtick.minor.width' : 2,
          'xtick.major.pad'   : 10,
          'xtick.major.width' : 4,
          'xtick.major.size'  : 12,
	  'ytick.minor.size'  : 6,
	  'ytick.minor.width' : 2,
          'ytick.major.width' : 4,
          'ytick.major.size'  : 12,
          'xtick.labelsize'    : 23,
          'ytick.labelsize'    : 23,
          'text.usetex'        : True,
	  'font.size'          : 30,
	  'font'               : 'serif',
          'mathtext.default'  : 'serif:bold'}
pl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts}"]
pl.rcParams.update(params)
pl.rc("font",family="serif")


colormap = ['b','g','r','m','c','y','DarkOrange', 'DarkBlue', 'Chartreuse']

kernen = [0]

for i in kernen:

     filename = "./relative/{:s}_tb_3d.txt".format(nuc[i])
     X = np.loadtxt(filename,unpack=False)
     
     fig = pl.figure(figsize=(6,6))
     fig.subplots_adjust(bottom=0.22,left=0.22)
     ax  = fig.add_subplot(111)
     
     
     k  = X[:,0] 
     P = X[:,1]
     MF_X = X[:,2]
     
     
     MF_norm = 2.0/(A[i]*(A[i]-1)) # because in convention corr normed on 1, MF will be normed lower than 1, correct for this, divide by A[i] to get normed to A again
     

     grid_x, grid_y = np.mgrid[0:max(k):800j, 0:max(P):800j]

     values = interpolate.griddata((k, P), MF_X, (grid_x, grid_y), method='linear')
     

     ax.imshow(values, origin='lower', extent=[0, max(P), 0, max(k)],
               aspect='auto', vmin=1e-4, vmax=10)
     pl.colorbar()
     
     
     ax.set_xlabel(r"$P_{12}$ [fm$^{-1}$]")
     ax.set_ylabel(r"$k_{12}$ [fm$^{-1}$]")
     
     pl.savefig("{:s}_tb_3d.png".format(nuc[i]), dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)

