import numpy as np
import pylab as pl
import scipy.integrate
###############################
### global variables   ########
###############################

nuc = ["Be","O","Al","Ar","Ca","Fe","Ag","Pb"]
A   = ["{9}","{16}","{27}","{40}","{40}","{56}","{109}","{208}"]
Z   = ["{4}","{8}","{13}","{18}","{20}","{26}","{47}","{82}"]

NORMALISATION_CONVENTION = "1" # choose "A" or "1", yes those are strings!



#################################
# set a bunch of plot parameters#
#################################

params = {'legend.fontsize' : 25,
          'legend.linewidth': 2,
          'axes.linewidth'  : 3,
          'axes.labelsize'  : 26,
	  'xtick.minor.size'  : 3,
	  'xtick.minor.width' : 0.2,
          'xtick.major.pad'   : 6,
          'xtick.major.width' : 4,
          'xtick.major.size'  : 9,
	  'ytick.minor.size'  : 5,
	  'ytick.minor.width' : 1.2,
          'ytick.major.width' : 3,
          'ytick.major.size'  : 8,
          'xtick.labelsize'    : 17,
          'ytick.labelsize'    : 17,
          'text.usetex'        : True,
	  'font.size'          : 30,
	  'font'               : 'serif',
          'mathtext.default'  : 'serif:bold'}
pl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{tensor} \usepackage{amsfonts}"]
pl.rcParams.update(params)
pl.rc("font",family="serif")


colormap = ['b','g','r','m','c','y','DarkOrange', 'DarkBlue', 'Chartreuse']
fig = pl.figure(figsize=(13,17))
fig.subplots_adjust(bottom=0.22,left=0.22, 
                    wspace = 0.08,  
                    hspace = 0.08)
ax_tot = fig.add_subplot(111)
ax_tot.axis('off')

for i in range(0,len(nuc)-2):

     filename = "{:s}_tb_rel.txt".format(nuc[i])
     X = np.loadtxt(filename,unpack=False)
     
     
     ax  = fig.add_subplot(3,2,i+1)
     ax.set_yscale('log')
     
     k  = X[0:200,0] 
     MF = X[0:200,1]
     if NORMALISATION_CONVENTION=="A":
          MF_norm = 1.0/A[i] # because in convention corr normed on 1, MF will be normed lower than 1, correct for this, divide by A[i] to get normed to A again
     elif NORMALISATION_CONVENTION=="1":
          MF_norm = 1 # because in convention corr normed on 1, MF will be normed lower than 1, correct for this
     else:
          raise Exception("Unsupported normalisation convention: \"{:s}\"".format(NORMALISATION_CONVENTION))
     ax.set_ylim((1e-4,1e3))
     ax.set_xlim((0,3))

     ax.plot(k,MF*MF_norm,color='b',label=r"{:s}".format(nuc[i]),lw=3)         
     ax.legend(frameon=False,ncol=3,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
     
     if((i+1) % 2 != 0):
          ax.set_ylabel(r"$n^{[2]}(k_{12})$ [fm$^{3}$]")
     else:
          ax.set_yticklabels([])
     if(i == len(nuc)-4 or i == len(nuc)-3):
          ax.set_xlabel(r"$k_{12}$ [fm$^{-1}$]")
     else:
          ax.set_xticklabels([])

pl.savefig("multi_tb_rel.png", dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
