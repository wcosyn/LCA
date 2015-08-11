import numpy as np
import pylab as pl
import scipy.integrate
###############################
### global variables   ########
###############################

nuc = ["He","Be","C","Al","Ar","Ca40","Fe","Ag","Pb"]
A   = [4,9,12,27,40,40,56,109,208]
NORMALISATION_CONVENTION = "A" # choose "A" or "1", yes those are strings!


#################################
# set a bunch of plot parameters#
#################################

params = {'legend.fontsize' : 25,
          'legend.linewidth': 2,
          'axes.linewidth'  : 3.5,
          'axes.labelsize'  : 40,
	  'xtick.minor.size'  : 6,
	  'xtick.minor.width' : 2,
          'xtick.major.pad'   : 10,
          'xtick.major.width' : 4,
          'xtick.major.size'  : 12,
	  'ytick.minor.size'  : 6,
	  'ytick.minor.width' : 2,
          'ytick.major.width' : 4,
          'ytick.major.size'  : 12,
          'xtick.labelsize'    : 33,
          'ytick.labelsize'    : 33,
          'text.usetex'        : True,
	  'font.size'          : 30,
	  'font'               : 'serif',
          'mathtext.default'  : 'serif:bold'}
pl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts}"]
pl.rcParams.update(params)
pl.rc("font",family="serif")


colormap = ['b','g','r','m','c','y','DarkOrange', 'DarkBlue', 'Chartreuse']

for i in range(len(nuc)):

     filename = "{:s}_ob_mf.txt".format(nuc[i])
     X = np.loadtxt(filename,unpack=False)
     fig = pl.figure(figsize=(6,6))
     fig.subplots_adjust(bottom=0.22,left=0.22)
     ax  = fig.add_subplot(111)
     ax.set_yscale('log')
     
     k  = X[0:29,0] # cut from [0:25] because a lot of noise for MF for high k values!
     MF = X[0:29,1]
     if NORMALISATION_CONVENTION=="A":
          MF_norm = 1.0/A[i] # because in convention corr normed on 1, MF will be normed lower than 1, correct for this, divide by A[i] to get normed to A again
     elif NORMALISATION_CONVENTION=="1":
          MF_norm = scipy.integrate.trapz(MF*k**2.,k) # because in convention corr normed on 1, MF will be normed lower than 1, correct for this
     else:
          raise Exception("Unsupported normalisation convention: \"{:s}\"".format(NORMALISATION_CONVENTION))
     
     ax.plot(k,MF/MF_norm,color='b',label=nuc[i],lw=3)         
     ax.legend(frameon=False,ncol=3,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
     ax.set_xlabel(r"$p$ [fm$^{-1}$]")
     ax.set_ylabel(r"$n^{[1]}(p)$ [fm$^{3}$]")
     ax.set_ylim((1e-4,1e3))
     pl.savefig("{:s}_ob_mf.png".format(nuc[i]), dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)

