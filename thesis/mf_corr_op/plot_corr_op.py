import numpy as np
import pylab as pl
import scipy.integrate
###############################
### global variables   ########
###############################

nuc = ["He","Be","C","O","Al","Ca40","Ca48","Fe","Ag"]
A   = [4,9,12,16,27,40,48,56,109]
NORMALISATION_CONVENTION = "1" # choose "A" or "1", yes those are strings!
OUTPUT_MF=True # True for output, False if you don't want output files

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

fig = pl.figure()
fig.subplots_adjust(bottom=0.22,left=0.22)
ax  = fig.add_subplot(111)
ax.set_yscale('log')

for i in range(len(nuc)):

     filename = "./{:s}/dens_ob4.00.111.{:s}.-1-1-1-1".format(nuc[i],nuc[i])
     X = np.loadtxt(filename,unpack=False)

     ax.plot(X[:,0],X[:,3],label=nuc[i],lw=3)
     k  = X[0:25,0] # cut from [0:25] because a lot of noise for MF for high k values!
     MF = X[0:25,1]
     if NORMALISATION_CONVENTION=="A":
          MF_norm = scipy.integrate.trapz(MF*k**2.,k)/A[i] # because in convention corr normed on 1, MF will be normed lower than 1, correct for this, divide by A[i] to get normed to A again
     elif NORMALISATION_CONVENTION=="1":
          MF_norm = scipy.integrate.trapz(MF*k**2.,k) # because in convention corr normed on 1, MF will be normed lower than 1, correct for this
     else:
          raise Exception("Unsupported normalisation convention: \"{:s}\"".format(NORMALISATION_CONVENTION))
     
     f = open("{:s}_mf_paper.txt".format(nuc[i]), "w")
     for j in range(len(k)):
          f.write(str(k[j]) + '\t' + str(MF[j]/MF_norm) + '\n')
     f.close()
          

     ax.plot(k,MF/MF_norm,'k--',lw=3)         
     print("# MF   normalisation of nucleus {:s} = {:e} ".format(nuc[i], scipy.integrate.trapz(MF/MF_norm*k**2.,k)))
     print("# Corr normalisation of nucleus {:s} = {:e} ".format(nuc[i], scipy.integrate.trapz(X[:,3]*X[:,0]**2,X[:,0])))
     if OUTPUT_MF:
          np.savetxt("MF_{:s}_normed.dat".format(nuc[i]),np.transpose(np.concatenate(([k],[MF/MF_norm]),axis=0)))

ax.legend(frameon=False,ncol=3,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
ax.set_xlabel(r"$p$ [fm$^{-1}$]")
ax.set_ylabel(r"$n^{[1]}(p)$ [fm$^{3}$]")
ax.set_ylim((1e-4,5e2))
pl.show()
