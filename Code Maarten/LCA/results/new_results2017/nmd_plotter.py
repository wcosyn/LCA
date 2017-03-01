#
# this is a script to plot/compare the results for 
# the momentum distributions I get from Maartens code
# 
# If you specify two files, it will compare those files
#
# use as python nmd_plotter.py file1 [file2]
#
# input files you want are probably named similar to
# 
# dens_ob4.00.111.He4.-1-1-1-1
# dens_rel.00.111.O.-1-1-1-1.-1-1
# 
# for the naming conventions: there is probably a README.txt
# file in the directory the data files are in explaining it.
#

import numpy as np
import pylab as pl
import sys
import os

#################################
# set a bunch of plot parameters#
#################################

params = {'legend.fontsize' : 16,
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




def ob_makeplot(filename1,nucleus):
    X1 = np.loadtxt(filename1,unpack=False)

    fig = pl.figure()
    fig.subplots_adjust(left=0.2,bottom=0.2)
    fig.suptitle("{:s}".format(nucleus))
    ax = fig.add_subplot(111)
     
    ax.plot(X1[:,0],X1[:,3],color='black',lw=2,label="mf+corr")
    ax.plot(X1[:,0],X1[:,1],color='green',lw=2,label="mf")
    ax.plot(X1[:,0],X1[:,2],color='firebrick',lw=2,label="corr")
    ax.plot(X1[:,0],X1[:,4],color='Goldenrod',lw=2,ls='dashed',label="central")
    ax.plot(X1[:,0],X1[:,5],color='blue',lw=2,ls='dashed',dashes=[1,3,2,3],label="tensor")
    ax.plot(X1[:,0],X1[:,6],color='teal',lw=2,ls='dashed',dashes=[4,2,6,3],label="spin/iso")

    ax.set_xlabel(r"$\mathbf{p \;\; [\text{fm}^{-1}}]$")
    ax.set_ylabel(r"$\mathbf{n^{[1]}(p) \;\; [\text{fm}^{3}}]$")
    ax.legend(frameon=False,numpoints=1,labelspacing=0)
    ax.set_yscale('log')
    ax.set_ylim((1e-5,ax.get_ylim()[1]))
    
    pl.savefig("obnmd_{:s}.pdf".format(nucleus))
    pl.show()

def tb_rel_makeplot(filename1,nucleus):
    X1 = np.loadtxt(filename1,unpack=False)

    fig = pl.figure()
    fig.subplots_adjust(left=0.2,bottom=0.2)
    fig.suptitle("{:s}".format(nucleus))
    ax = fig.add_subplot(111)
     
    ax.plot(X1[:,0],X1[:,4],color='black',lw=2,label="mf+corr")
    ax.plot(X1[:,0],X1[:,1],color='green',lw=2,label="mf")
    ax.plot(X1[:,0],X1[:,2],color='firebrick',lw=2,label="corr")
    ax.plot(X1[:,0],X1[:,3],color='chocolate',lw=2,label="3b corr")
    ax.plot(X1[:,0],X1[:,5],color='Goldenrod',lw=2,ls='dashed',label="central")
    ax.plot(X1[:,0],X1[:,6],color='blue',lw=2,ls='dashed',dashes=[1,3,2,3],label="tensor")
    ax.plot(X1[:,0],X1[:,7],color='teal',lw=2,ls='dashed',dashes=[4,2,6,3],label="spin/iso")

    ax.set_xlabel(r"$\mathbf{k_{12} \;\; [\text{fm}^{-1}}]$")
    ax.set_ylabel(r"$\mathbf{n^{[2]}(k_{12}) \;\; [\text{fm}^{3}}]$")
    ax.legend(frameon=False,numpoints=1,labelspacing=0)
    ax.set_yscale('log')
    ax.set_ylim((1e-5,ax.get_ylim()[1]))
    
    pl.savefig("tbnmd_rel_{:s}.pdf".format(nucleus))
    pl.show()
    

def ob_compareplot(filename1,filename2,nucleus):

    X1 = np.loadtxt(filename1,unpack=False)
    X2 = np.loadtxt(filename2,unpack=False)


    fig = pl.figure()
    fig.subplots_adjust(left=0.2,bottom=0.2)
    fig.suptitle("{:s}".format(nucleus))
    ax = fig.add_subplot(111)
    
    ax.plot([0],[0],mec='black',marker='x',mew=1,ms=8,lw=0,label="new calcs.")
    
    ax.plot(X2[:,0],X2[:,4],mec='black',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,4],color='black',lw=2,label="mf+corr")
    
    ax.plot(X2[:,0],X2[:,1],mec='green',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,1],color='green',lw=2,label="mf")

    ax.plot(X2[:,0],X2[:,2],mec='firebrick',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,2],color='firebrick',lw=2,label="corr")
    
    ax.plot(X2[:,0],X2[:,3],mec='chocolate',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,3],color='chocolate',lw=2,label="3b corr")

    ax.plot(X2[:,0],X2[:,5],mec='Goldenrod',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,5],color='Goldenrod',lw=2,ls='dashed',label="central")

    ax.plot(X2[:,0],X2[:,6],mec='blue',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,6],color='blue',lw=2,ls='dashed',dashes=[1,3,2,3],label="tensor")
    
    ax.plot(X2[:,0],X2[:,7],mec='teal',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,7],color='teal',lw=2,ls='dashed',dashes=[4,2,6,3],label="spin/iso")

    ax.set_xlabel(r"$\mathbf{k_{12} \;\; [\text{fm}^{-1}}]$")
    ax.set_ylabel(r"$\mathbf{n^{[2]}(k_{12}) \;\; [\text{fm}^{3}]}$")
    ax.legend(frameon=False,numpoints=1,labelspacing=0,loc=3)
    ax.set_yscale('log')
    ax.set_ylim((1e-11,ax.get_ylim()[1]))
    
    pl.savefig("obnmd_comparison_{:s}.pdf".format(nucleus))

def tb_rel_compareplot(filename1,filename2,nucleus):

    X1 = np.loadtxt(filename1,unpack=False)
    X2 = np.loadtxt(filename2,unpack=False)


    fig = pl.figure()
    fig.subplots_adjust(left=0.2,bottom=0.2)
    fig.suptitle("{:s}".format(nucleus))
    ax = fig.add_subplot(111)
    
    ax.plot([0],[0],mec='black',marker='x',mew=1,ms=8,lw=0,label="new calcs.")
    
    ax.plot(X2[:,0],X2[:,3],mec='black',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,3],color='black',lw=2,label="mf+corr")
    
    ax.plot(X2[:,0],X2[:,1],mec='green',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,1],color='green',lw=2,label="mf")

    ax.plot(X2[:,0],X2[:,2],mec='firebrick',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,2],color='firebrick',lw=2,label="corr")

    ax.plot(X2[:,0],X2[:,4],mec='Goldenrod',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,4],color='Goldenrod',lw=2,ls='dashed',label="central")

    ax.plot(X2[:,0],X2[:,5],mec='blue',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,5],color='blue',lw=2,ls='dashed',dashes=[1,3,2,3],label="tensor")
    
    ax.plot(X2[:,0],X2[:,6],mec='teal',marker='x',mew=1,ms=8,lw=0)
    ax.plot(X1[:,0],X1[:,6],color='teal',lw=2,ls='dashed',dashes=[4,2,6,3],label="spin/iso")

    ax.set_xlabel(r"$\mathbf{p \;\; [\text{fm}^{-1}}]$")
    ax.set_ylabel(r"$\mathbf{n^{[1]}(p) \;\; [\text{fm}^{3}}]$")
    ax.legend(frameon=False,numpoints=1,labelspacing=0,loc=3)
    ax.set_yscale('log')
    ax.set_ylim((1e-11,ax.get_ylim()[1]))
    
    pl.savefig("tbnmd_rel_comparison_{:s}.pdf".format(nucleus))
    pl.show()

def ob_compareResults():
    newresultfolder = "new_results2017/obnmd/"
    oldresultfolder = "result_1408"
    for filename in filter(lambda x: x.startswith("dens_ob4"),os.listdir(newresultfolder)):
        nucleus  = filename.split(".")[3] # in the naming convention of Maarten nucleus name will be between 3rd and 4th dot
        print("Nucleus is {:s}".format(nucleus))
        if os.path.isfile(os.path.join(oldresultfolder,nucleus,filename)):
            ob_compareplot(os.path.join(newresultfolder,filename),os.path.join(oldresultfolder,nucleus,filename),nucleus)
        else:
            ob_makeplot(filename)

def getNucleusFromFilename(f):
    nuc = os.path.basename(f).split(".")[3].replace("_",".") # in the naming convention of Maarten nucleus name will be between 3rd and 4th dot, replace underscores because tex won't like them
    print("[info] filename is {:s} , extracted nucleus is {:s}".format(f,nuc))
    return nuc 

if __name__=="__main__":
    compareplot = None
    makeplot    = None
    if ("ob" in sys.argv[1]):
        print("[info] ob detected from filename")
        compareplot = ob_compareplot
        makeplot    = ob_makeplot
    if ("rel" in sys.argv[1]):
        print("[info] tb rel detected from filename")
        compareplot = tb_rel_compareplot
        makeplot    = tb_rel_makeplot

    if len(sys.argv) == 2: # argument supplied, filename of data
        filename = sys.argv[1]
        makeplot(filename,getNucleusFromFilename(filename))
    elif len(sys.argv) == 3: # 2 files supplied, compare them
        print("#[Info] solid lines are from the first file you supplied")
        f1 = sys.argv[1]
        f2 = sys.argv[2]
        compareplot(f1,f2,getNucleusFromFilename(f1))
    else:
        compareResults()
