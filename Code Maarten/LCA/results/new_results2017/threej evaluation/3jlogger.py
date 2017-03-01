import numpy as np
import pylab as pl
import sys
import re
import os

#
# this is a script to investigate
# the 3j-coefficients usage
# give filename as parameter
# to process.
# searches for lines tagged with
# [3jlogger]
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
          'xtick.minor.size'  : 4,
          'ytick.minor.width' : 2,
          'ytick.minor.size'  : 4,
          'xtick.labelsize'    : 28,
          'ytick.labelsize'    : 28,
          'text.usetex'        : True,
          'font.size'          : 30}
#          'font'               : 'serif',
#          'mathtext.default'  : 'serif:bold'}
pl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts}"]
pl.rcParams.update(params)
pl.rc("font",family="serif")


def parse3jline(line):
    s3j = re.search('\(.*?\)',line).group(0)
    n =  int(re.search('::\s*(.*)\s*:',line).group(1))
    return s3j,n

def getNucFromFilename(fname):
    A = int(re.search(r'\d+',fname).group(0))
    nuc = re.search(r't(\D+)\d',fname).group(1)
    return A,nuc


def parsefile(fname):
    f = open(fname)
    A,nuc = getNucFromFilename(fname)
    m3j = {}
    ns = []
    maxj = 0
    maxj3j = None
    for line in f:
        if line.startswith("[3jlogger]"):
            t3j,n = parse3jline(line)
            assert( t3j not in m3j )
            tjm = max(map(int,t3j.strip("()").split(",")))
            if (tjm > maxj):
                maxj = tjm
                maxj3j = t3j
            m3j[t3j] = n
            ns.append(n)
    print("[Info][{:s}] total number of different 3js {:d}".format(nuc,len(maxj3j)))
    print("[Info][{:s}] total number of 3j evals is {:d}".format(nuc,sum(ns)))
    print("[Info][{:s}] maximum j value is  {:d} --> {:s}".format(nuc,maxj,maxj3j))
    return A,nuc,m3j,ns,maxj,maxj3j


def makesingle():
    A,nuc,m3j,ns,maxj,maxj3j = parsefile(sys.argv[1])
    nuclabel = r"$^{{ {:d} }}${:s}".format(A,nuc)
    fig = pl.figure(figsize=(12,7))
    fig.subplots_adjust(bottom=0.2,left=0.15)
    ax = fig.add_subplot(111)
    ax.hist(ns,bins=np.logspace(-0.1,int(np.log10(max(ns)))+1,100),rwidth=0.8,facecolor='lightseagreen')
    #ax.hist(ns,bins=np.linspace(0.5,max(ns)+0.5,max(ns)+1),rwidth=0.5,facecolor='lightseagreen')
    #ax.set_yscale('log')
    ax.set_xlim((0.5,ax.get_xlim()[1]))
    ax.set_xscale('log')
    print("maximum of ns is {:d}".format(max(ns)))
    print("number of 1 evals is {:d}".format( ns.count(1)))
    ax.set_ylim((0.5,ax.get_ylim()[1]))
    ax.set_yscale('log',nonposy='clip')
    ax.set_ylabel(r"\# of 3j coefficients ")
    ax.set_xlabel(r"\# evaluations (ob)")
    ax.annotate(nuclabel+"\n"+r"$2j_{{\text{{max}}}}={:d}$".format(maxj),xy=(0.75,0.75),xycoords='axes fraction')

def makestacked():
    files = filter( lambda n: re.match("t.*\.txt",n),os.listdir('.'))
    files.sort(key = lambda f: getNucFromFilename(f)[0])
    fig = pl.figure(figsize=(10,5))
    fig.subplots_adjust(bottom=0.2,left=0.15)
    ax = fig.add_subplot(111)
    nn = []
    nl = []
    As = []
    for fname in files:
        A,nuc,m3j,ns,maxj,maxj3j = parsefile(fname)
        nuclabel = r"$^{{ {:d} }}${:s}".format(A,nuc)
        nn.append(ns)
        nl.append(nuclabel)
        As.append(A)
    ax.hist(nn,bins=np.linspace(0.5,max(ns)+0.5,max(ns)+1),rwidth=0.8,stacked=False,label=nl,fill=True)
    ax.legend(frameon=False,fontsize=16)
    ax.set_yscale('log',nonposy='clip')
    ax.set_ylabel(r"\# of 3j coefficients ")
    ax.set_xlabel(r"\# evaluations (ob)")

    fig = pl.figure()
    fig.subplots_adjust(left=0.16,bottom=0.18)
    ax = fig.add_subplot(111)
    ax.plot(As, [ sum(ni) for ni in nn],'ro')
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(r"$A$")
    ax.set_ylabel(r"total \# of 3j evaluations")


if len(sys.argv) > 1:
    makesingle()
else:
    makestacked()
pl.show()

