#!/usr/bin/python
import sys
import numpy as np
import pylab as pl
import matplotlib.cm as cm
from matplotlib import rc
rc('text', usetex=True)
rc('font', size=24)
#rc('font',family='serif')
nuc = ["He","Be","C","O","Al","Ca40","Ca48","Fe","Ag","Ar","Pb"]
A   = [4,9,12,16,27,40,48,56,109,40,208]

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

len= len(sys.argv)-1
fig=pl.figure()
fig.subplots_adjust(bottom=0.15, left=0.15)
filename = sys.argv[2]
nuclide = filename[0] + filename[1]
if(nuclide[1] == "_"):
     nuclide = filename[0]
fig.text(.5,.75,nuclide,
        horizontalalignment='center',
        transform=fig.transFigure)


data = np.loadtxt(sys.argv[2], usecols=[0,1] )
X = data[:,0]
previous = data[:,1]
pl.fill_between( X*np.sqrt(2), 0, previous, color=cm.hot(1./(len) ))
pl.plot(X*np.sqrt(2), previous, color=cm.hot(1./(len)), label=  r"$l=0$" )
pl.xlim((0,4))

for i in range( 3, len ):
  print i
  print sys.argv[i]
  new = np.loadtxt( sys.argv[i], usecols=[1] )
  next = previous+ new
  pl.fill_between( X*np.sqrt(2), previous, next, color=cm.hot(float(i-1)/(len) ), label="l=%d" % (i-2) )
  pl.plot(X*np.sqrt(2), next, color=cm.hot( float(i-1)/(len) ), label=r"$l\leq%d$" % (i-2) )
  previous = next

total = np.loadtxt( sys.argv[len], usecols=[1] )
pl.plot( X*np.sqrt(2), total , color="black", label=r"total")
pl.legend()
#pl.xlabel(r" $P =  \frac{1}{\sqrt{2}} ( p_1+p_2 )$")
pl.xlabel(r" $P =  ( p_1+p_2 )\ [fm^{-1}]$")
pl.ylabel( r"$n^{[2]}(P) \ [fm^{3}]$" )
normalization = sum( total*X*X) * (X[1]-X[0] )
print "Normalization is ", normalization

#pl.show()
pl.savefig(sys.argv[1], dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)



