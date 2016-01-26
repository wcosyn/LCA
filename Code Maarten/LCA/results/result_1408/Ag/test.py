# A test file for checking if I can make sense of Maartens results...
import numpy as np
import pylab as pl

# columns in these files are k,mf,corr,total,cen,ten,siso,ce/te,ce/si,te/si
# where ce is central, te is tensor and si is spin isospin (I guess??)
Xnn = np.loadtxt("dens_ob4.-1-1.111.Ag.-1-1-1-1")
Xnp = np.loadtxt("dens_ob4.-11.111.Ag.-1-1-1-1")
Xpp = np.loadtxt("dens_ob4.11.111.Ag.-1-1-1-1")

XNN = np.loadtxt("dens_ob4.00.111.Ag.-1-1-1-1")


fig = pl.figure()
ax1 = fig.add_subplot(121)
ax1.plot(Xnn[:,0],Xnn[:,1],'g-')
ax1.plot(Xnn[:,0],Xnn[:,2],'g-',lw=3,label='nn')

ax1.plot(Xnp[:,0],Xnp[:,1],'r-')
ax1.plot(Xnp[:,0],Xnp[:,2],'r-',lw=3,label='np')

ax1.plot(Xpp[:,0],Xpp[:,1],'b-')
ax1.plot(Xpp[:,0],Xpp[:,2],'b-',lw=3,label='pp')

ax1.set_yscale('log')
ax1.set_ylim((1e-10,200))
ax1.legend()

ax2 = fig.add_subplot(122)
ax2.plot(XNN[:,0],XNN[:,1],'k-')
ax2.plot(XNN[:,0],XNN[:,2],'k-',lw=3,label='00')

ax2.plot(Xnn[:,0],Xnn[:,1]+Xnp[:,1]+Xpp[:,1],'r-')
ax2.plot(Xnn[:,0],Xnn[:,2]+Xnp[:,2]+Xpp[:,2],'r-',lw=3,label='nn+np+pp')

ax2.set_yscale('log')
ax2.set_ylim((1e-10,200))
ax2.legend()


from scipy.integrate import trapz

MF_norm_nn = trapz(Xnn[:,1]*Xnn[:,0]**2,Xnn[:,0])
MF_norm_np = trapz(Xnp[:,1]*Xnp[:,0]**2,Xnp[:,0])
MF_norm_pp = trapz(Xpp[:,1]*Xpp[:,0]**2,Xpp[:,0])

print("MF nn/np is {:}, expected N(N-1)/(2NZ) {:}".format(MF_norm_nn/MF_norm_np, 61.*60./(2.*61.*47.)))
print("MF pp/np is {:}, expected N(N-1)/(2NZ) {:}".format(MF_norm_pp/MF_norm_np, 47.*46./(2.*61.*47.)))

pl.show()
