import numpy as np
import pylab as pl
import collections
import matplotlib.patches
import scipy.optimize
import periodictable
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
          'xtick.labelsize'    : 28,
          'ytick.labelsize'    : 28,
          'text.usetex'        : True,
          'font.size'          : 30}
#          'font'               : 'serif',
#          'mathtext.default'  : 'serif:bold'}
pl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts}"]
pl.rcParams.update(params)
pl.rc("font",family="serif")


A,Z,t = np.loadtxt("./timeinfo.txt",unpack=True)
N=A-Z
shells = [2,6,8,14,16,20,28,32,38,40,50,58,64,68,70,82,92,100,106,110,112,126]
# shellnaming convention nlj, l \in {s,p,d,f,etc} 
# l = 0 <-> s
# l = 1 <-> p
# l = 2 <-> d
# l = 3 <-> f
# l = 4 <-> g
# l = 5 <-> h
# l = 6 <-> i
shellnames = ["0s1/2","0p3/2","0p1/2","0d5/2",
              "1s1/2","0d3/2","0f7/2","1p3/2",
              "0f5/2","1p1/2","0g9/2","0g7/2",
              "1d5/2","1d3/2","2s1/2","0h11/2",
              "0h9/2","1f7/2","1f5/2","2p3/2",
              "2p1/2","0i13/2"]

assert( len(shellnames) == len(shells))

Zvalshell = np.searchsorted(shells,Z)
Nvalshell = np.searchsorted(shells,N)
numshells = Zvalshell + Nvalshell + 2 # + 2 because start counting from zero
shelldict = collections.defaultdict(list) # group nuclei with same valence shells in dict with valence shells as key
for i in range(len(A)):
    zshell = Zvalshell[i]
    nshell = Nvalshell[i]
    print("Z,N = ({:3d},{:3d}) -- val. shells = ({:2d},{:2d}) = ({:5s},{:5s})".format(int(Z[i]),int(N[i]),zshell,nshell,shellnames[zshell],shellnames[nshell]))
    shelldict[(zshell,nshell)].append(i)



fig = pl.figure()
fig.subplots_adjust(bottom=0.2,left=0.16)
ax  = fig.add_subplot(111)
fexp = lambda x,a,b : a*np.exp(b*x)
fpow = lambda x,a,b : a*x**b
popt1,pcov1 = scipy.optimize.curve_fit(fexp,A,t,p0=(4,0.1))
popt2,pcov2 = scipy.optimize.curve_fit(fpow,A,t)

for shells,indices in shelldict.items():
    ts = map( lambda i: t[i], indices)
    As = map( lambda i: A[i], indices)
    Zs = map( lambda i: Z[i], indices)
    padding = 1.3 # log padding! (division and multiplication)
    x,y = min(As)/padding,min(ts)/padding
    w,h = max(As)*padding - x, max(ts)*padding - y
    ax.add_patch( matplotlib.patches.Rectangle( (x,y),w,h,fill=False) )
    shnames = map( lambda i: shellnames[i], shells)
    ax.annotate('({:},{:})'.format(*shnames),xy=(x,y+h),xytext=(0.45*x/padding,y*padding*2.3),fontsize=10,
                arrowprops=dict(facecolor='black',width=0.1,frac=0.,headwidth=0))
    elements = map( lambda i: periodictable.elements[Z[i]][A[i]],indices)
    elementnames = ",".join( map( lambda el: r"$^{{ {:d} }}${:s}".format(el.isotope,el.element), elements) )
    ax.annotate(elementnames,xy=(x+w,y),xytext=(2.3*x*padding,y/padding*0.45),fontsize=10,
                arrowprops=dict(facecolor='black',width=0.1,frac=0.,headwidth=0))


Ar = np.linspace(4,208,100)
ax.plot(A,t,marker='o',color='white',mec='black',mew=1,ms=8,lw=0)
ax.plot(Ar,fpow(Ar,*popt2),ls='dashed',lw=1,color='green')
#ax.plot(Ar,fexp(Ar,*popt1),ls='dashed',lw=3,color='red')
print("Extrapolation to A=208 using powerlaw fit = {:.2e} s = {:.2f} days".format(fpow(208.,*popt2),fpow(208.,*popt2)/86400.))
print("Power of A = {:.2f}".format(popt2[1]))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r"$A$")
ax.set_ylabel(r"comp. time (s)")


fig = pl.figure()
fig.subplots_adjust(bottom=0.2,left=0.16)
ax = fig.add_subplot(111)
ax.plot(Zvalshell,t,marker='o',color='red',mec='black',mew=1,ms=8,lw=0)
ax.plot(Nvalshell,t,marker='o',color='white',mec='black',mew=1,ms=8,lw=0)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r"number of shells")
ax.set_ylabel(r"comp. time (s)")

fig = pl.figure()
fig.subplots_adjust(bottom=0.2,left=0.16)
ax = fig.add_subplot(111)
ax.plot(numshells,t,marker='o',lw=0)

f = lambda x,a,b: a*x**b
popt,pcov = scipy.optimize.curve_fit(f,numshells,t)
ttemp = np.linspace(4,38,38-4+1)
ax.plot(ttemp,f(ttemp,*popt),'k--')
print("Pb 208 has 16+22=38 shells -> estimated computer time is {:.2e} s = {:.2f} days ".format(f(38,*popt),f(38,*popt)/86400.))
print("power of extrapolation : {:.2f}".format(popt[1]))


ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel("number of shells")
ax.set_ylabel(r"comp. time (s)")

pl.show()
