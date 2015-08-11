import numpy as np
import math
from distribution import Distribution
from scipy import special
import pylab as pl
import scipy.stats

#############################
## parameters for plotting ##
#############################


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



#########################
## nuclide parameters ##
########################

nuc = ["He","Be","C","O","Al","Ca40","Ca48","Fe","Ag","Ar","Pb"]
A   = [4,9,12,16,27,40,48,56,109,40,208]
Z   = [2,4,6,8,13,20,20,26,47,18,82]
N = [m - n for m,n in zip(A,Z)] #number of neutrons

Omega = []
for i in range(len(A)):
     Omega.append((45.0*(A[i]**(-1./3.)))-(25.0*(A[i]**(-2./3.))))
     
hbarc = 197.327 #(MeV*fm)
M_n = 939.565 #mass neutron (MeV) [Particle Data Group]
M_p = 938.272 #mass proton (MeV) [Particle Data Group]
M_average = (M_n + M_p)/2.



#############################
## radial HO wave function ##
#############################
  
#returns HO constant \nu for different mass
def v(M, index):
     return M*Omega[index] / (hbarc**2)  
     
#returns HO constant \nu in momentum space for different mass
def v_mom(M, index):
     return (hbarc**2) / (Omega[index]*M) 

#returns unnormalized radial wavefunction for HO calculated with constant a 
def Rad_wavefunc(k, a, n, l):
	t = a*(k**2)
	return (k**l)*special.eval_genlaguerre(n, l+0.5, t)*np.exp(-(a/2.)*(k**2))
     
#norm radial wave function
def norm(a, n, l):
	return np.sqrt((2*math.factorial(n)*(a**(l+3./2.)))/math.gamma(n+l+(3./2.)))


#orbital probability distribution for one-patricle momentum k
def prob(k):
     n = 0 
     l = 1 
     return 0.5*((Rad_wavefunc(k, v_mom(M_average,2), n, l)*norm(v_mom(M_average,2), n, l)*k)**2 + (Rad_wavefunc(k, v_mom(M_average,2), n, 0)*norm(v_mom(M_average,2), n, 0)*k)**2)



################
## simulation ##
################

#number of draws from distribution
size = 2000000


#theta and phi: draw from uniform distribution on a sphere
u1 = np.random.uniform(0.0, 1.0, size)
v1 = np.random.uniform(0.0, 1.0, size)

u2 = np.random.uniform(0.0, 1.0, size)
v2 = np.random.uniform(0.0, 1.0, size)

th1 = np.arccos(2*u1-1)
ph1 = v1*2*np.pi

th2 = np.arccos(2*u2-1)
ph2 = v2*2*np.pi

#one-particle momentum distribution
xr  = np.linspace(0,6,size+1)
pdf = prob(xr)
dist = Distribution(pdf,transform=Distribution.xmap(xr))

k1 = np.array(dist(size))
k2 = np.array(dist(size)) 

veck1 = np.array([k1*np.cos(ph1)*np.sin(th1),k1*np.sin(ph1)*np.sin(th1),k1*np.cos(th1)])
veck2 = np.array([k2*np.cos(ph2)*np.sin(th2),k2*np.sin(ph2)*np.sin(th2),k2*np.cos(th2)])


#relative 2body momentum vector and magnitude
vecK = (1./(2**0.5))*(veck1-veck2)
K = np.sum(np.multiply(vecK,vecK),axis=0)**0.5

#CoM 2body momentum vector and magnitude
vecP = (1./(2**0.5))*(veck1+veck2)
P = np.sum(np.multiply(vecP,vecP),axis=0)**0.5

#angle between relative momentum and CoM  momentum
theta_KP = np.arccos(np.sum(vecK*vecP,axis=0)/(K*P))

###################
## make a figure ##
###################


fig = pl.figure(figsize=(16,8))

ax2 = fig.add_subplot(121)


ax2.set_xlabel(r"$|\vec{k}_1|\ [fm^{-1}]$")
ax2.set_ylabel(r"$n^{[1]}(|\vec{k}_1|)|\vec{k}_1|^2\ [fm] $")
ax2.set_xlim((0,3))
ax2.set_ylim((0,1.5))
ax2.text(2.3,1.2,"C")

ax2.hist(k1,normed=True,bins=32)
ax2.plot(xr,prob(xr),color='green',lw=3)

ax = fig.add_subplot(122)


data, edges = np.histogram(theta_KP,normed=True,bins=np.linspace(0.0,np.pi,32))
av = []
for i in xrange(0,len(edges)-1):
     av.append((edges[i] + edges[i+1])/2.)

data = data/np.sin(av) 

data_th = np.loadtxt("C_angle.txt",unpack=False)

ax.bar(av-((edges[1]-edges[0])/2),data,edges[1]-edges[0])
ax.plot(data_th[:,0],data_th[:,1],color='green',lw=3)
ax.text(2.5,0.55,"C")

#ax.legend(frameon=False,ncol=3,handlelength=1.2,labelspacing=0.2,handletextpad=0.1,columnspacing=0.0)
ax.set_xlabel(r"$\theta_{kP}$")
ax.set_ylabel(r"$n^{[2]}(\theta_{kP})$")
ax.set_xlim((0,np.pi))
ax.set_ylim((0,0.7))



pl.savefig("carbon",dpi=None,facecolor='w', edgecolor='w',orientation='portrait',papertype=None,transparent=False,bbox_inches='tight',pad_inches=0.1,frameon=None)
#pl.show()





