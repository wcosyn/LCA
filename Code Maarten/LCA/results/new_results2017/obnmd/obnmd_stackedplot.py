#
# make a plot with obnmd for different 
# nuclei
#
# launch with -suppressmf to suppress mf above k=3.5
#
# python obnmd_stackedplot.py [-suppressmf]


import numpy as np
import pylab as pl
import sys
import re
import os

#
# this is a script to make a plot
# of the one-body momentum distributions
# for a variety of nuclei
#

#################################
# set a bunch of plot parameters#
#################################

params = {'legend.fontsize' : 30,
#          'legend.linewidth': 2,
          'axes.linewidth'  : 3.5,
          'axes.labelsize'  : 30,
          'xtick.major.pad' : 10,
          'xtick.major.width' : 3,
          'xtick.major.size'  : 6,
          'ytick.major.width' : 3,
          'ytick.major.size'  : 6,
          'xtick.minor.width' : 2,
          'xtick.minor.size'  : 3,
          'ytick.minor.width' : 2,
          'ytick.minor.size'  : 3,
          'xtick.labelsize'    : 20,
          'ytick.labelsize'    : 20,
          'text.usetex'        : True,
          'font.size'          : 30}
#          'font'               : 'serif',
#          'mathtext.default'  : 'serif:bold'}
pl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts}"]
pl.rcParams.update(params)
pl.rc("font",family="serif")


#
# make a list of the nuclei you want to plot
#
nucs = ["He4","C12","Al27","Ar40","Ca40","Ca48","Fe56","Kr84_8cores_nondiscint","Ag108","Xe124","Nd142nondiscint","W184_16cores","Au197nondiscint","Pb208"]
filenames = [ "dens_ob4.00.111.{:s}.-1-1-1-1s".format(n) for n in nucs ]
nucl = map( lambda s: re.search(r"^[A-Z,a-z]+",s).group(0), nucs)
A    = map( lambda s: int(re.search(r"\d+",s).group(0)), nucs)
nuclabels = [ r"$^{{ \mathbf{{ {:d} }} }}${:s}".format(A[i],nucl[i]) for i in range(len(nucs)) ]

fig = pl.figure()
fig.subplots_adjust(bottom=0.15,left=0.15)
ax  = fig.add_subplot(111)
figk4 = pl.figure()
figk4.subplots_adjust(bottom=0.15,left=0.15)
axk4  = figk4.add_subplot(111)
for i in range(len(nucs)):
    print("#[info] reading file {:s}".format(filenames[i]))
    X = np.loadtxt(filenames[i],unpack=False)
    print("#[info] max rel diff of tot and mf+corr {:.2e}".format(np.max( np.abs( (X[:,1] + X[:,2] - X[:,3]) / X[:,3]))))
    if "-suppressmf" in sys.argv:
        kmax = 3.5
        print("#[info] suppressing mf above {:}".format(kmax))
        ki = np.searchsorted(X[:,0],kmax)
        X[ki:,1] = 0.
    ax.plot(X[:,0],X[:,1]+X[:,2],label=nuclabels[i],lw=3)
    axk4.plot(X[:,0],(X[:,1]+X[:,2])*X[:,0]**4,label=nuclabels[i],lw=3)

ax.set_yscale("log")
ax.legend(frameon=False,fontsize=16,ncol=2)

axk4.set_yscale("log")
axk4.legend(frameon=False,fontsize=16,ncol=2,loc=4)
ax.set_xlabel(r"$\mathbf{k}$ [fm$^{\mathbf{-1}}$]")
ax.set_ylabel(r"$\mathbf{ n^{[1]}(k)} $ [fm$^{\mathbf{3}}$]")

axk4.set_xlabel(r"$\mathbf{k}$ [fm$^{\mathbf{-1}}$]")
axk4.set_ylabel(r"$\mathbf{ k^{4} \;  n^{[1]}(k)} $ [fm$^{\mathbf{3}}$]")

fig.savefig("obnmd_stackedplot.pdf")
figk4.savefig("obnmd_stackedplot_k4.pdf")

pl.show()
