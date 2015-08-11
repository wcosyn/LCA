import numpy as np
import pylab as pl
import scipy.integrate
###############################
### global variables   ########
###############################

nuc = ["He","C","O","Al","Ar","Ca40","Fe","Ag","Pb"]
A   = [4,12,16,27,40,40,56,109,208]
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

     filename_X = "./results/{:s}_tb_rel.txt".format(nuc[i])
     X = np.loadtxt(filename_X,unpack=False)
     filename_Y = "./results_maarten/rel_dens_mf.{:s}".format(nuc[i])
     Y = np.loadtxt(filename_Y,unpack=False)
     
     #fig = pl.figure(figsize=(6,6))
     #fig.subplots_adjust(bottom=0.22,left=0.22)
     #ax  = fig.add_subplot(111)
     
     
     k_X  = X[0:150,0] 
     MF_X = X[0:150,1]
     MF_Y = Y[0:150,1]
     dist_Y = MF_Y*k_X*k_X
     norm_Y = scipy.integrate.trapz(dist_Y,x=None,dx = 0.02)
     print MF_X/(MF_Y/norm_Y)
      
     """ax.plot(k_X,MF_X,color='b',label=nuc[i],lw=3)
     ax.plot(k_X,MF_Y/norm_Y,color='black',linestyle='--',label=nuc[i],lw=3)
     ax.legend(frameon=False,ncol=3,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
     ax.set_xlabel(r"$k_{12}$ [fm$^{-1}$]")
     ax.set_ylabel("n^{rel}_2")
     ax.set_ylim((0.6,1.4))
     pl.savefig("./figures/{:s}_tb_mf.png".format(nuc[i]), dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)"""

