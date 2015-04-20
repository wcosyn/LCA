import numpy as np
import pylab as pl
import scipy.integrate
###############################
### global variables   ########
###############################



NORMALISATION_CONVENTION = "A" # choose "A" or "1", yes those are strings!


#################################
# set a bunch of plot parameters#
#################################

params = {'legend.fontsize' : 20,
          'legend.linewidth': 2,
          'axes.linewidth'  : 3.5,
          'axes.labelsize'  : 32,
	  'xtick.minor.size'  : 6,
	  'xtick.minor.width' : 2,
          'xtick.major.pad'   : 10,
          'xtick.major.width' : 4,
          'xtick.major.size'  : 12,
	  'ytick.minor.size'  : 6,
	  'ytick.minor.width' : 2,
          'ytick.major.width' : 4,
          'ytick.major.size'  : 12,
          'xtick.labelsize'    : 28,
          'ytick.labelsize'    : 28,
          'text.usetex'        : True,
	  'font.size'          : 30,
	  'font'               : 'serif',
          'mathtext.default'  : 'serif:bold'}
pl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts}"]
pl.rcParams.update(params)
pl.rc("font",family="serif")


colormap = ['b','g','r','m','c','y','DarkOrange', 'DarkBlue', 'Chartreuse']


filename_X = "wave_10.txt"
data_X = np.loadtxt(filename_X,unpack=False)
X = data_X[:,0]
M_X = data_X[:,1]

filename_Y = "wave_02.txt"
data_Y = np.loadtxt(filename_Y,unpack=False)
Y = data_Y[:,0]
M_Y = data_Y[:,1]

fig = pl.figure(figsize=(6,6))
fig.subplots_adjust(bottom=None,left=None)
ax  = fig.add_subplot(111)



ax.plot(X,M_X,color='b',lw=3, label="n = 1, l = 0") 
ax.plot(Y,M_Y,color='g',lw=3,label="n = 0, l = 2")         
ax.set_xlabel(r"$r[fm]$")
ax.set_ylabel(r"$R^2(r)r^2$")
ax.set_xlim((0,9))
ax.legend(frameon=False,ncol=1,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
pl.savefig("waves.png", dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)