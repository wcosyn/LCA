#!/usr/bin/python
import numpy as np
import pylab as pl
import os

int= 0
X= np.loadtxt( "cen/dens_ob3.11.100.Ag11.0.p", usecols=(0,) )
p_M= np.zeros_like(X)
p_C= np.zeros_like(X)
p_A= np.zeros_like(X)
n_M= np.zeros_like(X)
n_C= np.zeros_like(X)
n_A= np.zeros_like(X)

for E in range( 0, 9 ):
  new_M= np.zeros_like(X)
  new_C= np.zeros_like(X)
  new_A= np.zeros_like(X)

  for set1 in range( 1, 4, 1 ):
    for set2 in range( set1, 4, 1 ):
      file= "cen/dens_ob3.11.100.Ag{:d}{:d}.{:d}.p".format( set1, set2, E)
      if not os.path.isfile(file):
        print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
      else:
        X,M,C,A = np.loadtxt(file, unpack=True)
        new_M += M
        new_C += C
        new_A += A
        p_M += M
        p_C += C
        p_A += A

  test_int= (X[1]-X[0])* sum( X* X* new_A )
  int += test_int

  fname= "dens_ob3.11.100.Ag.{:d}.p".format(E)
  print fname
  new_file= open( fname, "w" )
  for i in range( 0, len(X) ):
    new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], new_M[i], new_C[i], new_A[i] ) )
  new_file.close()







  new_M= np.zeros_like(X)
  new_C= np.zeros_like(X)
  new_A= np.zeros_like(X)

  for set1 in range( 1, 4, 1 ):
    for set2 in range( set1, 5, 1 ):
      file= "cen/dens_ob3.-11.100.Ag{:d}{:d}.{:d}.p".format( set1, set2, E)
      if not os.path.isfile(file):
        print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
      else:
        X,M,C,A = np.loadtxt(file, unpack=True)
        new_M += M
        new_C += C
        new_A += A
        p_M += M
        p_C += C
        p_A += A

  test_int= (X[1]-X[0])* sum( X* X* new_A )
  int += test_int

  fname= "dens_ob3.-11.100.Ag.{:d}.p".format(E)
  print fname
  new_file= open( fname, "w" )
  for i in range( 0, len(X) ):
    new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], new_M[i], new_C[i], new_A[i] ) )
  new_file.close()

  new_M= np.zeros_like(X)
  new_C= np.zeros_like(X)
  new_A= np.zeros_like(X)

  for set1 in range( 1, 4, 1 ):
    for set2 in range( set1, 5, 1 ):
      file= "cen/dens_ob3.-11.100.Ag{:d}{:d}.{:d}.n".format( set1, set2, E)
      if not os.path.isfile(file):
        print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
      else:
        X,M,C,A = np.loadtxt(file, unpack=True)
        new_M += M
        new_C += C
        new_A += A
        n_M += M
        n_C += C
        n_A += A

  test_int= (X[1]-X[0])* sum( X* X* new_A )
  int += test_int

  fname= "dens_ob3.-11.100.Ag.{:d}.n".format(E)
  print fname
  new_file= open( fname, "w" )
  for i in range( 0, len(X) ):
    new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], new_M[i], new_C[i], new_A[i] ) )
  new_file.close()

  new_M= np.zeros_like(X)
  new_C= np.zeros_like(X)
  new_A= np.zeros_like(X)

  for set1 in range( 1, 5, 1 ):
    for set2 in range( set1, 5, 1 ):
      file= "cen/dens_ob3.-1-1.100.Ag{:d}{:d}.{:d}.n".format( set1, set2, E)
      if not os.path.isfile(file):
        print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
      else:
        X,M,C,A = np.loadtxt(file, unpack=True)
        new_M += M
        new_C += C
        new_A += A
        n_M += M
        n_C += C
        n_A += A

  test_int= (X[1]-X[0])* sum( X* X* new_A )
  int += test_int

  fname= "dens_ob3.-1-1.100.Ag.{:d}.n".format(E)
  print fname
  new_file= open( fname, "w" )
  for i in range( 0, len(X) ):
    new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], new_M[i], new_C[i], new_A[i] ) )
  new_file.close()



print int

fname= "dens_ob3.00.100.Ag.-1.n"
print fname
new_file= open( fname, "w" )
for i in range( 0, len(X) ):
  new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], n_M[i], n_C[i], n_A[i] ) )
new_file.close()

fname= "dens_ob3.00.100.Ag.-1.p"
print fname
new_file= open( fname, "w" )
for i in range( 0, len(X) ):
  new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], p_M[i], p_C[i], p_A[i] ) )
new_file.close()

fname= "dens_ob3.00.100.Ag.-1"
print fname
new_file= open( fname, "w" )
for i in range( 0, len(X) ):
  new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], p_M[i]+n_M[i], p_C[i]+n_C[i], p_A[i]+n_A[i] ) )
new_file.close()




