import numpy as np
import pylab as pl
import sys
import re

#
# this is a script to investigate
# the 3j-coefficients hash map performance
# give filename as parameter
# to process that file.
# searches for lines tagged with
# [3jmap]
#

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
          'xtick.minor.width' : 2,
          'xtick.minor.size'  : 3,
          'ytick.minor.width' : 2,
          'ytick.minor.size'  : 3,
          'xtick.labelsize'    : 20,
          'ytick.labelsize'    : 20,
          'text.usetex'        : True,
          'font.size'          : 30}
#          'font'               : 'serif',
#          'mathtext.default'  : 'serif:bold'}
pl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts}"]
pl.rcParams.update(params)
pl.rc("font",family="serif")

def parse3jline(line):
    n =  int(re.search('has\s*(\d+)\s*elements',line).group(1))
    return n

def parsefile(fname):
    f = open(fname)
    ns = []
    for line in f:
        if line.startswith("[3jmap] bucket"):
            ns.append(parse3jline(line))
    return ns


def getNucFromFilename(fname):
    A = int(re.search(r'\d+',fname).group(0))
    nuc = re.search(r't(\D+)\d',fname).group(1)
    return A,nuc

def main():
    ns = parsefile(sys.argv[1]) # ns is the number of elements in each bucket
    A,nuc = getNucFromFilename(sys.argv[1])
    nuclabel = r"$^{{ {:d} }}${:s}".format(A,nuc)
    el = sum(ns) # el is the number of elements in the map
    Nb = max(ns) #Nb is the maximum number of elements in a bucket
    fig = pl.figure()
    fig.subplots_adjust(left=0.18,bottom=0.2)
    ax = fig.add_subplot(111)
    bins = np.linspace(-0.5,Nb+0.5,Nb+2)
    n,_,_ = ax.hist(ns,bins=bins,rwidth=0.65,facecolor='firebrick',lw=2,ec='black')
    ax.set_xticks( np.linspace(0,Nb,Nb+1))
    ax.set_yscale('log',nonposy='clip')
    ax.set_ylabel(r"\# of buckets")
    ax.set_xlabel(r"bucket size")
    ax.annotate(nuclabel,xy=(0.75,0.75),xycoords='axes fraction')
    n = list(map(int,n)) # make n ints
    ncoll = sum(list(filter(lambda x: x>1,ns))) # count number of collisions

    totb = sum(n)
    zerob = n[0]
    oneb  = n[1]
    multb = sum(n[2:])
    print("[info] total number of elements : {:8d}".format(el))
    print("[info] total number of collisons: {:8d}  = {:.4f}".format(ncoll,1.*ncoll/el))
    print("[info] total number of buckets  : {:8d}".format(totb))
    print("[info] empty buckets            : {:8d}  = {:.4f}, :: random exp = {:.4f}".format(zerob,1.*zerob/totb,(1.-1./totb)**(el)))
    print("[info] single el. buckets       : {:8d}  = {:.4f}".format(oneb,1.*oneb/totb))
    print("[info] mult. el. buckets        : {:8d}  = {:.4f}".format(multb,1.*multb/totb))
    for i in range(2,len(n)):
        print("[info]  |---> buckets with {:d} el. : {:8d} = {:.4f}".format(i,n[i],1.*n[i]/totb))

    pl.show()

main()
