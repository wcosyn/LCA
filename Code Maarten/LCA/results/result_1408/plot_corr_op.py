import numpy as np
import pylab as pl
import scipy.integrate

nuc = ["He","Be","C","O","Al","Ca40","Ca48","Fe","Ag"]
A   = [4,9,12,16,27,40,48,56,109]

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
pl.rcParams.update(params)

fig = pl.figure()
fig.subplots_adjust(bottom=0.22,left=0.22)
ax  = fig.add_subplot(111)
ax.set_yscale('log')
for i in range(len(nuc)):
        # in the dens_ob4... files, columns are k, mf, corr, total, central , tensor, spinisospin, ce/te, ce/si, te/si
	X = np.loadtxt("./{:s}/dens_ob4.00.111.{:s}.-1-1-1-1".format(nuc[i],nuc[i]),unpack=False)
        normMf   = scipy.integrate.simps(X[:,1]*X[:,0]**2,X[:,0]) # remember, norm is times k**2
        normCorr = scipy.integrate.simps(X[:,2]*X[:,0]**2,X[:,0]) # remember, norm is times k**2
        normTot  = scipy.integrate.simps(X[:,3]*X[:,0]**2,X[:,0]) # remember, norm is times k**2
        print("#[Info] normalisation for {:4s},A={:4d}: (mf,corr,tot) = ({:.2e},{:.2e},{:.2e}), mf+corr-tot={: .2e}".format(nuc[i],A[i],normMf,normCorr,normTot,normMf+normCorr-normTot))
	ax.plot(X[:,0],X[:,3],label=nuc[i],lw=3)#/A[i])
	ax.plot(X[0:25,0],X[0:25,1],'k--',lw=3)#/A[i])
ax.legend(frameon=False,ncol=3,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
ax.set_xlabel(r"$p$ [fm$^{-1}$]")
ax.set_ylabel(r"$n^{[1]}(p)$ [fm$^{3}$]")
ax.set_ylim((1e-4,5e2))
pl.show()
