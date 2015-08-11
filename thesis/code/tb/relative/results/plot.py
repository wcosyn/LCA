#!/usr/bin/python
import sys
import numpy as np
import pylab as pl
import matplotlib.cm as cm
from matplotlib import rc
from scipy import integrate
rc('text', usetex=True)
rc('font', size=24)
#rc('font',family='serif')
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


len= len(sys.argv)-1
fig=pl.figure()
fig.subplots_adjust(bottom=0.22,left=0.22)
ax  = fig.add_subplot(111)

data = np.loadtxt( sys.argv[2], usecols=[0,1] )
X = data[:,0]
MF_X = data[:,1]
normalization = integrate.trapz(MF_X*X*X,X)
print "Normalization is ", normalization
ax.plot(X, MF_X, color='black' ,lw=3, label= sys.argv[3])
ax.set_yscale('log')

ax.legend(frameon=False,ncol=1,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
ax.set_xlabel(r"$k_{12}$ [fm$^{-1}$]")
ax.set_ylabel(r"$n^{[2]}(k_{12})$ [fm$^{3}$]")
ax.set_xlim((0,3))
ax.set_ylim((1e-4,5e1))
     


#pl.show()
pl.savefig(sys.argv[1], dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)



