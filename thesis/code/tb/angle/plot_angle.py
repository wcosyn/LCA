import numpy as np
import pylab as pl
import scipy.integrate
import sys

nuclide = sys.argv[1]
###############################
### global variables   ########
###############################

nuc = ["He","Be","C","O","Al","Ar","Ca40","Fe","Ag","Pb"]
A   = [4,9,12,16,27,40,40,56,109,208]
NORMALISATION_CONVENTION = "A" # choose "A" or "1", yes those are strings!

NUC = 0
if nuclide in nuc:
	for i in xrange(0,len(nuc)):
		if(nuclide == nuc[i]):
			NUC = i
			break
else:
	print "nuclide not found"

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



filename_X = "./new_results/{:s}_angle.txt".format(nuc[NUC])
X = np.loadtxt(filename_X,unpack=False)
     
     
fig = pl.figure(figsize=(6,6))
fig.subplots_adjust(bottom=0.22,left=0.22)
ax  = fig.add_subplot(111)
    
     
k_X  = X[:,0]
MF_X = (X[:,1]-X[:,1].min())/(X[:,1].max()-X[:,1].min())
     
     
     
     
MF_norm = 2.0/(A[NUC]*(A[NUC]-1)) # because in convention corr normed on 1, MF will be normed lower than 1, correct for this, divide by A[i] to get normed to A again
     
     
ax.plot(k_X,MF_X,color='black',lw=3)
ax.yaxis.tick_left()
ax.xaxis.tick_bottom()
ax.text(2.7,0.84,nuclide,horizontalalignment='center')
ax.legend(frameon=False,ncol=3,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
ax.set_xlabel(r"$\theta_{kP}$")
ax.set_ylabel(r"$\frac{n^{[2]}-n^{[2]}_{min}}{n^{[2]}_{max}-n^{[2]}_{min}}$")
ax.set_xlim([0,np.pi])
ax.set_ylim([MF_X.min(),1.03])
pl.xticks([0, np.pi/4,np.pi/2, 3*np.pi/4, np.pi],
          [r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',r'$\frac{3\pi}{4}$', r'$\pi$'])
pl.savefig("./figures/{:s}_tb_angle.png".format(nuc[NUC]), dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)

