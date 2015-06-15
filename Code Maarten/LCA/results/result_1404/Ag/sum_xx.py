#!/usr/bin/python
import numpy as np
import pylab as pl
import os

int= 0
X= np.loadtxt( "dens_ob4.-1-1.110.Ag.0.n", usecols=(0,) )
pp_M= np.zeros_like(X)
pp_C= np.zeros_like(X)
pp_A= np.zeros_like(X)
nn_M= np.zeros_like(X)
nn_C= np.zeros_like(X)
nn_A= np.zeros_like(X)
pn_M= np.zeros_like(X)
pn_C= np.zeros_like(X)
pn_A= np.zeros_like(X)


for E in range( 0, 9 ):
###################################################################3
#
#       PP
#
##################################################################3
  ppE_M= np.zeros_like(X)
  ppE_C= np.zeros_like(X)
  ppE_A= np.zeros_like(X)

  file= "dens_ob4.11.110.Ag.{:d}.p".format(  E)
  file_cp= "dens_ob4.11.110.Ag.{:d}".format(  E)
  if not os.path.isfile(file):
    print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
    print "pp", E
  else:
    X,M,C,A = np.loadtxt(file, unpack=True, usecols=(0,1,2,3))
    pp_M += M
    pp_C += C
    pp_A += A
    ppE_M += M
    ppE_C += C
    ppE_A += A
    os.system( "cp -v " + file + " " + file_cp )
###################################################################3
#
#       PN
#
##################################################################3

  pnE_M= np.zeros_like(X)
  pnE_C= np.zeros_like(X)
  pnE_A= np.zeros_like(X)
  file= "dens_ob4.-11.110.Ag.{:d}.p".format(  E)
  if not os.path.isfile(file):
    print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
    print "pn", E
  else:
    X,M,C,A = np.loadtxt(file, unpack=True, usecols=(0,1,2,3))
    pn_M += M
    pn_C += C
    pn_A += A
    pnE_M += M
    pnE_C += C
    pnE_A += A


  file= "dens_ob4.-11.110.Ag.{:d}.n".format(  E)
  if not os.path.isfile(file):
    print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
    print "pn", E
  else:
    X,M,C,A = np.loadtxt(file, unpack=True, usecols=(0,1,2,3))
    pn_M += M
    pn_C += C
    pn_A += A
    pnE_M += M
    pnE_C += C
    pnE_A += A

  fname= "dens_ob4.-11.110.Ag.{:d}".format(  E)
  print fname
  new_file= open( fname, "w" )
  for i in range( 0, len(X) ):
    new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], pnE_M[i], pnE_C[i], pnE_A[i] ) )
  new_file.close()

###################################################################3
#
#       NN
#
##################################################################3
  nnE_M= np.zeros_like(X)
  nnE_C= np.zeros_like(X)
  nnE_A= np.zeros_like(X)

  file= "dens_ob4.-1-1.110.Ag.{:d}.n".format(  E)
  file_cp= "dens_ob4.-1-1.110.Ag.{:d}".format(  E)
  if not os.path.isfile(file):
    print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
    print "nn", E
  else:
    X,M,C,A = np.loadtxt(file, unpack=True, usecols=(0,1,2,3))
    nn_M += M
    nn_C += C
    nn_A += A
    nnE_M += M
    nnE_C += C
    nnE_A += A
    os.system( "cp -v " + file + " " + file_cp )

  fname= "dens_ob4.00.110.Ag.{:d}".format(  E)
  print fname
  new_file= open( fname, "w" )
  for i in range( 0, len(X) ):
    new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], 
      ppE_M[i]+ pnE_M[i]+ nnE_M[i], 
      ppE_C[i]+ pnE_C[i]+ nnE_C[i], 
      ppE_A[i]+ pnE_A[i]+ nnE_A[i] ) )
  new_file.close()



########################################3
#   OUTPUT 00
################################3

int= (X[1]-X[0])* sum( X* X* (pp_A+ pn_A+ nn_A ) )
print int

fname= "dens_ob4.11.110.Ag.-1"
print fname
new_file= open( fname, "w" )
for i in range( 0, len(X) ):
  new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], pp_M[i], pp_C[i], pp_A[i] ) )
new_file.close()

fname= "dens_ob4.-11.110.Ag.-1"
print fname
new_file= open( fname, "w" )
for i in range( 0, len(X) ):
  new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], pn_M[i], pn_C[i], pn_A[i] ) )
new_file.close()

fname= "dens_ob4.-1-1.110.Ag.-1"
print fname
new_file= open( fname, "w" )
for i in range( 0, len(X) ):
  new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], nn_M[i], nn_C[i], nn_A[i] ) )
new_file.close()


