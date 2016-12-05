# Trying to figure Maartens results out
# script written by Camille, not necessarily 
# correct interpretation of Maartens results :p


import numpy as np
import pylab as pl
import matplotlib.ticker
import scipy.integrate

#################################
# set a bunch of plot parameters#
#################################

params = {'legend.fontsize' : 30,
          'axes.linewidth'  : 3.5,
          'axes.labelsize'  : 30,
          'xtick.major.pad' : 10,
          'xtick.major.width' : 3,
          'xtick.major.size'  : 6,
          'ytick.major.width' : 3,
          'ytick.major.size'  : 6,
          'xtick.labelsize'    : 28,
          'ytick.labelsize'    : 28,
          'text.usetex'        : True,
          'font.size'          : 30}
pl.rcParams.update(params)



X    = np.loadtxt("C/dens_ob4.11.111.C.-1-1-1-1",unpack=False)
k    = X[:,0]
mf   = X[:,1]
corr = X[:,2]
tot  = X[:,3]

mfnorm   = scipy.integrate.simps(y=mf*k**2,x=k)
corrnorm = scipy.integrate.simps(y=corr*k**2,x=k)
totnorm  = scipy.integrate.simps(y=tot*k**2,x=k)

print("#[Info] mf   integral : {:.2e}".format(mfnorm))
print("#[Info] corr integral : {:.2e}".format(corrnorm))
print("#[Info] tot  integral : {:.2e}".format(totnorm))

print("#[Info] max diff between tot and corr+mf is {:.2e}".format(np.max(tot-corr-mf)))

corrDashes=[8,3,8,3]
totDashes =[2,3,3,2]


k*= 197. # times HBARC to go from fm^-1 to MeV/c
fig = pl.figure()
fig.subplots_adjust(left=0.18,bottom=0.2)
ax = fig.add_subplot(111)
ax.plot(k,1./mfnorm*mf,ls='dashed',dashes=[3,3,3,3],color='red',lw=3,label="mean field")
ax.plot(k,1./corrnorm*corr,ls='solid',lw=3,color='black',label="LCA")
#ax.plot(k,corr,ls='dashed',dashes=corrDashes,lw=3,label="corr")
#ax.plot(k,tot,ls='dashed',dashes=totDashes,lw=3,label="tot")
#ax.plot(k,tot/(corr+mf),ls='solid',lw=3,color='pink')

ax.set_yscale("log")
ax.set_ylim((1e-3,10))
ax.set_xlim((0,900))
ax.set_xlabel(r"$k \textrm{ [MeV/c]}$")
ax.set_ylabel(r"$n(k) \textrm{ [fm}^{-3}\textrm{]}$")
ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
ax.legend(frameon=False)

pl.show()
