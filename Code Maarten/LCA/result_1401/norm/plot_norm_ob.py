#!/usr/bin/python
from matplotlib import rc
import sys
import os
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.colors as plc

rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], 
                                  'monospace': ['Computer Modern Typewriter']})
rc('font', size=35 )
rc('text', usetex=True)
rc('axes', lw=3)
rc('xtick.major', size=10)
rc('ytick.major', size=10)
#rc('xtick.major', pad=30)
#rc('mathtext', default='regular')
rc('xtick', labelsize=35)
pl.rcParams['text.latex.preamble']=[r"\usepackage{bm} \usepackage{amsmath}"]

def save(path, ext='pdf', show=True, verbose=True):
    """Save a figure from pyplot.
 
    Parameters
    ----------
    path : string
        The path (and filename, without the extension) to save the
        figure to.
 
    ext : string (default='pdf')
        The file extension. This must be supported by the active
        matplotlib backend (see matplotlib.backends module).  Most
        backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.
 
    close : boolean (default=True)
        Whether to close the figure after saving.  If you want to save
        the figure multiple times (e.g., to multiple formats), you
        should NOT close it in between saves or you will have to
        re-plot it.
 
    verbose : boolean (default=True)
        Whether to print information about when and where the image
        has been saved.
 
    """
    
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
        directory = '.'
 
    # If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)
 
    # The final path to save to
    savepath = os.path.join(directory, filename)
 
    if verbose:
        print("Saving figure to '%s'..." % savepath),
 
    # Actually save the figure
    pl.savefig(savepath)
    
    # Close it
    if show:
     pl.show(block=False)
 
    if verbose:
        print("Done")

data= np.loadtxt( "norm_ob.110.%s.-1" % (sys.argv[1]), unpack=True )
refnn= data[5]
refnp= data[6]
refpp= data[7]
axis= data[1]

llabels = ['s','p','d','f','g','h']
axislist=[None] * (len(axis))

for i in range( 0, len(axis) ):
  int= axis[i]
  n= np.int_(np.floor(int/100))
  l= np.int_(np.floor((int-n*100)/10))
  j= np.int_((int-n*100-l*10))
  label = '$' + `n` + llabels[l] + '\\frac{' + `j` + '}{2}$'
  axislist[i]= label

xtics= [x+0.5 for x in range( len(axis) )]
ytics= [x+0.5 for x in range( len(axis) )]

size= np.sqrt( len(refnn) )
refnn= refnn.reshape( (size,size) )
refnp= refnp.reshape( (size,size) )
refpp= refpp.reshape( (size,size) )
refpn= refnp.transpose()
refpppn= np.concatenate( (refpp,refpn), axis=1 )
refnpnn = np.concatenate( (refnp,refnn), axis=1 )
ref = np.concatenate( (refpppn, refnpnn) )


bshow=False

f1= pl.figure()
f1.subplots_adjust(bottom=0.155, left=0.14)
ax1= f1.add_subplot(111)
norm=plc.Normalize(vmin=0,vmax=np.max(ref) )
pl.title("All E")
pl.pcolormesh(ref, norm=norm)
pl.colorbar()
pl.xticks( xtics, axislist )
pl.yticks( ytics, axislist )
pl.text( -0.85, size/2-0.1, "p", rotation="vertical" )
pl.text( -0.85, 3*size/2-0.1, "n", rotation="vertical" )
pl.text( size/2-0.1, -0.75, "p")
pl.text( 3*size/2-0.1, -0.75, "n")
save("all.%s" % sys.argv[1], show=bshow )

norm_ratio=plc.Normalize(vmin= 0, vmax=1, clip=False )
for i in range( 2, len(sys.argv) ):
  data= np.loadtxt( "norm_ob.110.%s.%s" % (sys.argv[1], sys.argv[i]), unpack=True )
  conn= data[5]
  conp= data[6]
  copp= data[7]
  conn= conn.reshape( (size,size) )
  conp= conp.reshape( (size,size) )
  copp= copp.reshape( (size,size) )
  copn= conp.transpose()
  copppn= np.concatenate( (copp,copn), axis=1 )
  conpnn = np.concatenate( (conp,conn), axis=1 )
  co = np.concatenate( (copppn, conpnn) )

  f2= pl.figure()
  f2.subplots_adjust(bottom=0.155, left=0.155)
  ax2= f2.add_subplot(111)
  pl.title("E=%s" %(sys.argv[i]))
  pl.pcolormesh(co, norm=norm)
  pl.colorbar()
  pl.xticks( xtics, axislist )
  pl.yticks( ytics, axislist )
  pl.text( -1, size/2-0.1, "p", rotation="vertical" )
  pl.text( -1, 3*size/2-0.1, "n", rotation="vertical" )
  pl.text( size/2-0.1, -0.75, "p")
  pl.text( 3*size/2-0.1, -0.75, "n")
  save("E.%s.%s" % (sys.argv[i] ,sys.argv[1]), show=bshow )

  f3= pl.figure()
  f3.subplots_adjust(bottom=0.155, left=0.155)
  ax3= f3.add_subplot(111)
  pl.title("Ratio E=%s" %(sys.argv[i]) )
  pl.pcolormesh(co/ref, norm=norm_ratio)
  pl.colorbar()
  pl.xticks( xtics, axislist )
  pl.yticks( ytics, axislist )
  pl.text( -1, size/2-0.1, "p", rotation="vertical" )
  pl.text( -1, 3*size/2-0.1, "n", rotation="vertical" )
  pl.text( size/2-0.1, -0.75, "p")
  pl.text( 3*size/2-0.1, -0.75, "n")
  save("ratio.%s.%s" % (sys.argv[i] ,sys.argv[1]), show=bshow )

if bshow:
  pl.show()




#mf=data[0]
#corr=data[1]
#print "mf", mf
#print "corr", corr
#size= len(corr)
#resize= np.sqrt(size)
#corr2= np.reshape( corr, (resize,resize) )
#print corr2

#norm=mpl.colors.Normalize(vmin= 0, vmax=1.5, clip=False )
#pl.pcolormesh(corr2, norm=norm)
#pl.colorbar()
#pl.show()
