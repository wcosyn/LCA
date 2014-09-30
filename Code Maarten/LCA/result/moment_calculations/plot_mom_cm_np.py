#!/usr/bin/python
#
# ./plot_mom_cm_np.py plot_mom_cm.eps c al fe pb 12 27 56 208
#
import sys
import numpy as np
import pylab as pl
import matplotlib.cm as cm
from matplotlib import rc


rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], 
                                  'monospace': ['Computer Modern Typewriter']})
rc('font', size=30)
rc('text', usetex=True)
rc('axes', lw=3)
rc('xtick.major', size=10)
rc('ytick.major', size=10)
rc('mathtext', default='regular')
rc('xtick', labelsize=30)

#pl.rcParams['axes.linewidth'] = 5 # set the value globally

fig=pl.figure(figsize=(12,12))
fig.subplots_adjust(bottom=0.07, left=0.1, right=0.975, top=0.98, wspace=0, hspace=0)
linestyles = [ '--', ':', '-.']

len= (len(sys.argv)-2)/2

n_plots = len
n_cols = 2
n_rows = (n_plots + 1) // n_cols
for plot_num in range(n_plots):
  print plot_num
  data0 = np.loadtxt( "hocommom_np_00.%s" % sys.argv[plot_num+2] , usecols=[0,1] )
  data1 = np.loadtxt( "hocommom_np_01.%s" % sys.argv[plot_num+2] , usecols=[0,1] )
  X = data0[:,0] *np.sqrt(2)
  previous = (data0[:,1]+data1[:,1])/np.sqrt(8)

  total0 = np.loadtxt( "hocommom_np_-10.%s" % sys.argv[plot_num+2] , usecols=[1] )/np.sqrt(8)
  total1 = np.loadtxt( "hocommom_np_-11.%s" % sys.argv[plot_num+2] , usecols=[1] )/np.sqrt(8)
  total= total0+total1
  normalization = sum( total*X*X) * (X[1]-X[0] )
  total = total/normalization
  previous = previous/normalization

  print "Normalization is ", normalization
  ax = pl.subplot(n_rows, n_cols, plot_num+1)
  if plot_num == 2 or plot_num == 3:
    ax.set_xticklabels([])
    pl.xlabel(r"$P_{12}$ [GeV]")
  else:
    ax.set_xticklabels([])
  if plot_num == 0 or plot_num == 2:
    ax.set_yticklabels(["",50,100,150,200,250,300])
    pl.ylabel( "$P_2(P_{12})$")
  else: 
    ax.set_yticklabels([])
  if plot_num==3:
    ax.set_xticklabels(["",0.1,0.2,0.3,0.4,0.5])
  if plot_num==2:
    ax.set_xticklabels([0.0,0.1,0.2,0.3,0.4,0.5])

  pl.fill_between( X, 0, previous, color='0.7', label=r"$l=0$" )
  pl.plot(X, previous, linestyle='-', lw=3, color='0.7', label=r"$l=0$" )
  pl.plot( X, total , ls='-', color="black", lw=3, label=r"all $l$")
#  pl.text(0.2,150, "%s" % sys.argv[plot_num+2])
  pl.xlim([0,0.5])
  pl.ylim([0.0,300])

  if(plot_num==0):
    legend=pl.legend(loc=1,title=r"$^{%s}$%s" % (sys.argv[plot_num+2+len], sys.argv[plot_num+2] ) )
  else:
    legend=pl.legend(loc=1,title=r"$^{%s}$%s" % (sys.argv[plot_num+2+len], sys.argv[plot_num+2] ) )
  

  for l in ax.get_xticklines() + ax.get_yticklines(): 
#    l.set_markersize(10) 
    l.set_markeredgewidth(3)
#  [i.set_linewidth(3) for i in ax.spines.itervalues()]




#pl.show()
pl.savefig( sys.argv[1] )


#for i in range( 3, len ):
#  print i
#  print sys.argv[i]
#  new = np.loadtxt( sys.argv[i], usecols=[1] )/normalization/np.sqrt(8)
#  next = previous+ new
#  #pl.fill_between( X*np.sqrt(2), previous, next, color=cm.hot(float(i-1)/(len) ), label="l=%d" % (i-2) )
#  pl.plot(X, next, linestyles[i-2], lw=2, color="black",label=r"$l\leq%d$" % (i-2) )
#  previous = next

##pl.xlabel(r" $P =  \frac{1}{\sqrt{2}} ( p_1+p_2 )$")





