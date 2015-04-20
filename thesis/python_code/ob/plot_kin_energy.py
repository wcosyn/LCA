import numpy as np
import pylab as pl
import scipy.integrate
###############################
### global variables   ########
###############################

nuc = ["O12","O13","O14","O15","O16","O17","O18","O19","O20","O21","O22","O23","O24"]
A   = [12,13,14,15,16,17,18,19,20,21,22,23,24]

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


filename = "kine_O_iso.txt"
X = np.loadtxt(filename,unpack=False, usecols = (2,))
fig = pl.figure(figsize=(6,6))
fig.subplots_adjust(bottom=0.22,left=0.22)
ax  = fig.add_subplot(111)

e = len(A)    
ke  = X

ax.plot(A,ke,color='b',lw=3)         
ax.set_xlabel(r"A")
ax.set_ylabel(r"$T_{kin}$ [MeV]")
pl.savefig("kin_energy_O_iso.png", dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)

