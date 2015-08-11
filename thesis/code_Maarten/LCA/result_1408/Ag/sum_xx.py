#!/usr/bin/python
import numpy as np
import pylab as pl
import os
import sys

A = sys.argv[1]

sel=""
if( len(sys.argv) > 2 ):
  sel= "."+ sys.argv[2]



int= 0
X= np.loadtxt( "dens_ob4.11.111.{:s}.-1-1-1-1{:s}".format(A,sel), usecols=(0,) )


###################################################################3
#
#       PP
#
##################################################################3

file= "dens_ob4.11.111.{:s}.-1-1-1-1{:s}".format(A,sel)
if not os.path.isfile(file):
  print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
  print "pp", E
else:
  Y_pp= np.loadtxt(file, unpack=False,)
###################################################################3
#
#       PN
#
##################################################################3
file= "dens_ob4.-11.111.{:s}.-1-1-1-1{:s}".format(A,sel)
if not os.path.isfile(file):
  print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
  print "pn", E
else:
  Y_pn= np.loadtxt(file, unpack=False,)

###################################################################3
#
#       NN
#
##################################################################3
file= "dens_ob4.-1-1.111.{:s}.-1-1-1-1{:s}".format(A,sel)
if not os.path.isfile(file):
  print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
  print "nn", E
else:
  Y_nn= np.loadtxt(file, unpack=False,)


Y = Y_pp+ Y_pn+ Y_nn

print Y

fname= "dens_ob4.00.111.{:s}.-1-1-1-1{:s}".format(A,sel)
print fname
new_file= open( fname, "w" )


for line in open("dens_ob4.11.111.{:s}.-1-1-1-1{:s}".format(A, sel)):
  if "norm" in line:
    new_file.write(line)
    norm= float(line.split("=")[-1])
    break



for i in range( 0, len(X) ):
  new_file.write("{:f}\t".format( X[i] ) )
  for j in range( 1, len(Y[i]) ):
    new_file.write("{:f}\t".format( Y[i][j] ) )
  new_file.write("\n")

new_file.close()



########################################3
#   OUTPUT 00
################################3

print Y.transpose()[3]
int= (X[1]-X[0])* sum( X* X* (Y.transpose())[3]  )
print int



