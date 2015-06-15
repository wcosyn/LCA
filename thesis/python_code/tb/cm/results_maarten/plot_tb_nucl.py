import numpy as np
import pylab as pl
import scipy.integrate
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

kernen = [2]

for i in kernen:

     filename_X = "{:s}_tot_maarten_normed.txt".format(nuc[i])
     X = np.loadtxt(filename_X,unpack=False)
     
     fig = pl.figure(figsize=(6,6))
     fig.subplots_adjust(bottom=0.22,left=0.22)
     ax  = fig.add_subplot(111)
     ax.set_yscale('log')
     
     k_X  = X[0:150,0] # cut from [0:25] because a lot of noise for MF for high k values!
     MF_X = X[0:150,1]
     
     
     MF_norm = 2.0/(A[i]*(A[i]-1)) # because in convention corr normed on 1, MF will be normed lower than 1, correct for this, divide by A[i] to get normed to A again
     
     
     ax.plot(k_X,MF_X,color='b',label=nuc[i],lw=3)
     
     ax.legend(frameon=False,ncol=3,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
     ax.set_xlabel(r"$P_{12} = k_1 + k_2$ [fm$^{-1}$]")
     ax.set_ylabel(r"$n^{[2]}(P_{12})$ [fm$^{3}$]")
     ax.set_ylim((1e-4,0.1e2))
     pl.savefig("{:s}_tb_cm_mf.png".format(nuc[i]), dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
