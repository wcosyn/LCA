#
# this is a script to compare the results for 
# the one-body momentum distributions I get from Maartens code
# against his old results.
#

import numpy as np
import pylab as pl
import sys
import os

#################################
# set a bunch of plot parameters#
#################################

params = {'legend.fontsize' : 20,
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





def compareplot(filename1,filename2,nucleus):

    X1 = np.loadtxt(filename1,unpack=False)
    X2 = np.loadtxt(filename2,unpack=False)


    fig = pl.figure()
    fig.subplots_adjust(left=0.2,bottom=0.2)
    ax = fig.add_subplot(111)
    ax.plot(X2[:,0],X2[:,1],color='white',marker='o',mew=2,mec='green',lw=0,label="{:s} mf new".format(nucleus))
    ax.plot(X1[:,0],X1[:,1],color='green',lw=2,label="{:s} mf old".format(nucleus))

    ax.plot(X2[:,0],X2[:,2],color='white',marker='o',mew=2,mec='firebrick',lw=0,label="{:s} corr new".format(nucleus))
    ax.plot(X1[:,0],X1[:,2],color='firebrick',lw=2,label="{:s} corr old".format(nucleus))

    ax.set_xlabel(r"$\mathbf{p \;\; \text{fm}^{-1}}$")
    ax.set_ylabel(r"$\mathbf{n^{[1]}(p) \;\; \text{fm}^{3}}$")
    ax.legend(frameon=False,numpoints=1)
    ax.set_yscale('log')
    ax.set_ylim((1e-5,ax.get_ylim()[1]))

    pl.savefig("obnmd_comparison_{:s}.pdf".format(nucleus))

    pl.show()

if __name__=="__main__":
    newresultfolder = "new_results2017/obnmd/"
    oldresultfolder = "result_1408"
    for filename in filter(lambda x: x.startswith("dens_ob4"),os.listdir(newresultfolder)):
        nucleus  = filename.split(".")[3] # in the naming convention of Maarten nucleus name will be between 3rd and 4th dot
        print("Nucleus is {:s}".format(nucleus))
        compareplot(os.path.join(newresultfolder,filename),os.path.join(oldresultfolder,nucleus,filename),nucleus)

