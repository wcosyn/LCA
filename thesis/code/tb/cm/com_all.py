import numpy as np
import pylab as pl
import scipy.integrate

HBARC = 197.327


A   = [4,12,56]

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

fig = pl.figure()
ax = fig.add_subplot(111)
for i in nuc : 
     
     k1,p1 = np.loadtxt("./results_maarten/" + i + "_cm_maarten_normed.txt",unpack=True)
     k2,p2 = np.loadtxt(i + "_tb_cm.txt",unpack=True)

     HBARC = 197.327
     k2 *= HBARC/1.e3

     n1 = scipy.integrate.simps(p1*k1**2,k1)
     n2 = scipy.integrate.simps(p2*k2**2,k2)

     print("Norm is   :  {:f}".format(n1))
     print("Norm is   :  {:f}".format(n2))
     print("Factor is :  {:f}".format(n1/n2))
     
     ax.plot(k2,p2*n1/n2,label=i,lw=3)
     ax.plot(k1,p1,color = 'black',linestyle = '--', label=i,lw=3)
     
     """k,p = np.loadtxt(i + "_tb_cm.txt",unpack=True)
     k *= HBARC/1.e3
     n = scipy.integrate.simps(p*k**2,k)
     ax.plot(k,p/n,label=i,lw=3)"""
     
     ax.legend(frameon=False,ncol=3,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
     ax.set_xlabel(r"$P_{12} = \frac{1}{\sqrt{2}} (k_1 + k_2)$ [$GeV$]")
     ax.set_ylabel(r"$n^{[2]}(P_{12})$ [$GeV^{-3}$]")
     ax.set_ylim((0,1050))
     ax.set_xlim((0,0.4))
pl.savefig("./figures/tb_cm_all.pdf", dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
     
