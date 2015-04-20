#!/usr/bin/python
import sys
import numpy as np
import pylab as pl
import matplotlib.cm as cm
from matplotlib import rc
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
previous = data[:,1]
ax.fill_between( X*np.sqrt(2), 0, previous, color=cm.hot(1./(len) ))
ax.plot(X*np.sqrt(2), previous, color=cm.hot(1./(len) ), label=r"$l=0$" )

for i in range( 3, len ):
  print i
  print sys.argv[i]
  new = np.loadtxt( sys.argv[i], usecols=[1] )
  next = previous+ new
  ax.fill_between( X*np.sqrt(2), previous, next, color=cm.hot(float(i-1)/(len) ), label="l=%d" % (i-2) )
  ax.plot(X*np.sqrt(2), next, color=cm.hot( float(i-1)/(len) ), label=r"$l\leq%d$" % (i-2) )
  previous = next

total = np.loadtxt( sys.argv[len], usecols=[1] )
ax.plot( X*np.sqrt(2), total , color="black", label=r"total")
ax.legend(frameon=False,ncol=1,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
#pl.xlabel(r" $P =  \frac{1}{\sqrt{2}} ( p_1+p_2 )$")
ax.set_xlabel(r"$P_{12} = k_1 + k_2$ [fm$^{-1}$]")
ax.set_ylabel(r"$n^{[2]}(P_{12})$ [fm$^{3}$]")
ax.set_xlim((0,4))
     
normalization = sum( total*X*X) * (X[1]-X[0] )
print "Normalization is ", normalization

#pl.show()
pl.savefig(sys.argv[1], dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)



