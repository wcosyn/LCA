#!/usr/bin/python
import sys
import numpy as np
import pylab as pl
import matplotlib.cm as cm
from matplotlib import rc

rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], 
                                  'monospace': ['Computer Modern Typewriter']})

## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

rc('font', size=35 )

rc('text', usetex=True)
rc('axes', lw=3)
rc('xtick.major', size=10)
rc('ytick.major', size=10)
rc('mathtext', default='regular')
rc('xtick', labelsize=35)
pl.rcParams['text.latex.preamble']=[r"\usepackage{bm} \usepackage{amsmath}"]
#pl.rcParams['axes.linewidth'] = 5 # set the value globally

fig=pl.figure(figsize=(12,12))
fig.subplots_adjust(bottom=0.07, left=0.1, right=0.975, top=0.98, wspace=0, hspace=0)
linestyles = [ '--', ':', '-.']

sys.argv= [ "", "plot_mom_cm_pp.pdf", "C", "Al", "Fe", "Pb", 12, 27, 56, 208]

size= (len(sys.argv)-2)/2

n_plots = size
n_cols = 2
n_rows = (n_plots + 1) // n_cols
for plot_num in range(n_plots):
  print plot_num
  data = np.loadtxt( "dens_com_pp.00-1.%s" % sys.argv[plot_num+2] , usecols=[0,1] )
  X = data[:,0]*0.197327
  previous = data[:,1]/0.197327/0.197327/0.197327

  data_total = np.loadtxt( "hocommom_pp_-1-1.%s" % sys.argv[plot_num+2] , usecols=[0,1] )
  X_total = data_total[:,0] *np.sqrt(2)
  total = data_total[:,1]/np.sqrt(8)
  normalization = sum( total*X_total*X_total) * (X_total[1]-X_total[0] )
  total = total/normalization
  previous = previous/normalization

  print "Normalization is ", normalization
  ax = pl.subplot(n_rows, n_cols, plot_num+1)
  if plot_num == 2 or plot_num == 3:
    ax.set_xticklabels([])
    pl.xlabel(r"$\bm{P_{12}\mathrm{[GeV]}}$")
  else:
    ax.set_xticklabels([])
  if plot_num == 0 or plot_num == 2:
    yticklabels= [' ',50,100,150,200,250,300]
    yticklabels_bold=[]
    for i in range( len(yticklabels) ):
      yticklabels_bold.append(r'$\bm{'+ str(yticklabels[i])+ '}$')
    ax.set_yticklabels(yticklabels_bold)
    pl.ylabel( r"$\bm{P_2(P_{12})}$")
  else: 
    ax.set_yticklabels([])
  if plot_num==3:
    xticklabels= ['',0.1,0.2,0.3,0.4,0.5]
    xticklabels_bold=[]
    for i in range( len(xticklabels) ):
      xticklabels_bold.append(r'$\bm{'+ str(xticklabels[i])+ '}$')
    ax.set_xticklabels(xticklabels_bold)
  if plot_num==2:
    xticklabels= [0.0,0.1,0.2,0.3,0.4,0.5]
    xticklabels_bold=[]
    for i in range( len(xticklabels) ):
      xticklabels_bold.append(r'$\bm{'+ `xticklabels[i]`+ '}$')
    ax.set_xticklabels(xticklabels_bold)

  pl.fill_between( X, 0, previous*5, color='0.7', label=r"$l=0$" )
  p1, = pl.plot(X, previous*5, linestyle='-', lw=3, color='0.7', label=r"$\displaystyle \mathrm{l=0}$" )
  p2, = pl.plot( X_total, total , ls='-', color="black", lw=3, label=r"all $\displaystyle \mathrm{l}$")
#  pl.text(0.2,150, "%s" % sys.argv[plot_num+2])
  pl.xlim([0,0.5])
  pl.ylim([0.0,300])

  abcd=['(a)','(b)','(c)','(d)']

#  if(plot_num==0):
  legend=pl.legend( [p1,p2], [ r"$\bm{nl}\bf{=}\bm{00\ \mathrm{[x5]}}$",r"$\bm{\mathrm{all}}\ \bm{nl}$"], bbox_to_anchor=[0.62, 0.99], frameon=False, loc=9,
      title=r"$\bm{\mathrm{%s}\; ^{%s}\mathrm{%s}}$" % (abcd[plot_num],sys.argv[plot_num+2+size], sys.argv[plot_num+2] ) )

# legend=pl.legend( [p1,p2], [ r"($nl$)=($00$) [x5]", "all ($nl$)"], bbox_to_anchor=[0.65, 0.99], frameon=False, loc=9,title=r"$^{%s}$%s" % (sys.argv[plot_num+2+size], sys.argv[plot_num+2] ) )
#  else:
#    legend=pl.legend([p1,p2], [ r"($nl$)=($00$)", "all ($nl$)"], bbox_to_anchor=[0.65, 0.95], loc=9,title=r"$^{%s}$%s" % (sys.argv[plot_num+2+size], sys.argv[plot_num+2] ) )
#    legend=pl.legend( [],[], frameon=False, bbox_to_anchor=[0.65, 0.99], loc=9,title=r"$^{%s}$%s" % (sys.argv[plot_num+2+size], sys.argv[plot_num+2] ) )
  

  for l in ax.get_xticklines() + ax.get_yticklines(): 
#    l.set_markersize(10) 
    l.set_markeredgewidth(3)
#  [i.set_linewidth(3) for i in ax.spines.itervalues()]




pl.savefig( sys.argv[1], bbox_inches='tight', pad_inches=0.05, format='pdf')
#pl.show()
print sys.argv[1]
#pl.savefig( sys.argv[1], bbox_inches='tight' )


#for i in range( 3, size ):
#  print i
#  print sys.argv[i]
#  new = np.loadtxt( sys.argv[i], usecols=[1] )/normalization/np.sqrt(8)
#  next = previous+ new
#  #pl.fill_between( X*np.sqrt(2), previous, next, color=cm.hot(float(i-1)/(size) ), label="l=%d" % (i-2) )
#  pl.plot(X, next, linestyles[i-2], lw=2, color="black",label=r"$l\leq%d$" % (i-2) )
#  previous = next

##pl.xlabel(r" $P =  \frac{1}{\sqrt{2}} ( p_1+p_2 )$")





