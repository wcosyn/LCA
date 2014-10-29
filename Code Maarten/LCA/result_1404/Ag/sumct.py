#!/usr/bin/python
import numpy as np
import pylab as pl
import os

int= 0
X= np.loadtxt( "raw_data/dens_ob4.-1-1.110.Ag11.0.n", usecols=(0,) )
p_T= np.zeros_like(X)
p_C= np.zeros_like(X)
p_CT= np.zeros_like(X)
n_T= np.zeros_like(X)
n_C= np.zeros_like(X)
n_CT= np.zeros_like(X)

for E in range( 0, 9 ):
  new_T= np.zeros_like(X)
  new_C= np.zeros_like(X)
  new_CT= np.zeros_like(X)

###################################################################3
#
#       PP
#
##################################################################3
  for set1 in range( 1, 5, 1 ):
    for set2 in range( set1, 5, 1 ):
      if set1 == 4 or set2 == 4:
        file= "raw_data/dens_ob4.11.110.Ag{:d}{:d}.{:d}.p".format( set1, set2, E)
      else:
        file= "raw_data/dens_ob4.-1-1.110.Ag{:d}{:d}.{:d}.n".format( set1, set2, E)

      if not os.path.isfile(file):
        print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
        print "pp", set1, set2
      else:
        X,C,T,CT = np.loadtxt(file, unpack=True, usecols=(0,4,5,7))
        new_T += T
        new_C += C
        new_CT += CT
        p_T += T
        p_C += C
        p_CT += CT

  test_int= (X[1]-X[0])* sum( X* X* new_C )
#  print test_int

  int += test_int

  fname= "dens_ob4.11.ct0.Ag.{:d}.p".format(E)
  print fname
  new_file= open( fname, "w" )
  for i in range( 0, len(X) ):
    new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], new_T[i], new_C[i], new_CT[i] ) )
  new_file.close()



###################################################################3
#
#       PN
#
##################################################################3




  new_T= np.zeros_like(X)
  new_C= np.zeros_like(X)
  new_CT= np.zeros_like(X)

  for set1 in range( 1, 5, 1 ):
    for set2 in range( set1, 7, 1 ):
      if set1 == 4 or set2 == 4:
        file= "raw_data/dens_ob4.-11.110.Ag{:d}{:d}.{:d}.p".format( set1, set2, E)
      else:
        file= "raw_data/dens_ob4.-11.110.Ag{:d}{:d}.{:d}.n".format( set1, set2, E)
      if not os.path.isfile(file):
        print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
      else:
        X,C,T,CT = np.loadtxt(file, unpack=True, usecols=(0,4,5,7))
        new_T += T
        new_C += C
        new_CT += CT
        p_T += T
        p_C += C
        p_CT += CT

  test_int= (X[1]-X[0])* sum( X* X* new_C )
  int += test_int

  fname= "dens_ob4.-11.ct0.Ag.{:d}.p".format(E)
  print fname
  new_file= open( fname, "w" )
  for i in range( 0, len(X) ):
    new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], new_T[i], new_C[i], new_CT[i] ) )
  new_file.close()

  new_T= np.zeros_like(X)
  new_C= np.zeros_like(X)
  new_CT= np.zeros_like(X)

  for set1 in range( 1, 5, 1 ):
    for set2 in range( set1, 7, 1 ):
      file= "raw_data/dens_ob4.-11.110.Ag{:d}{:d}.{:d}.n".format( set1, set2, E)
      if not os.path.isfile(file):
        print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
      else:
        X,C,T,CT=  np.loadtxt(file, unpack=True, usecols=(0,4,5,7))
        new_T += T
        new_C += C
        new_CT += CT
        n_T += T
        n_C += C
        n_CT += CT

  test_int= (X[1]-X[0])* sum( X* X* new_C )
#  print test_int

  int += test_int

  fname= "dens_ob4.-11.ct0.Ag.{:d}.n".format(E)
  print fname
  new_file= open( fname, "w" )
  for i in range( 0, len(X) ):
    new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], new_T[i], new_C[i], new_CT[i] ) )
  new_file.close()

  new_T= np.zeros_like(X)
  new_C= np.zeros_like(X)
  new_CT= np.zeros_like(X)

###################################################################3
#
#       NN
#
##################################################################3
  for set1 in range( 1, 7, 1 ):
    for set2 in range( set1, 7, 1 ):
      file= "raw_data/dens_ob4.-1-1.110.Ag{:d}{:d}.{:d}.n".format( set1, set2, E)
      if not os.path.isfile(file):
        print "ERROR FILE NOT FOUND, NOT POSSIBLE " , file
      else:
        X,C,T,CT = np.loadtxt(file, unpack=True, usecols=(0,4,5,7))
        new_T += T
        new_C += C
        new_CT += CT
        n_T += T
        n_C += C
        n_CT += CT

  test_int= (X[1]-X[0])* sum( X* X* new_C )
#  print test_int

  int += test_int

  fname= "dens_ob4.-1-1.ct0.Ag.{:d}.n".format(E)
  print fname
  new_file= open( fname, "w" )
  for i in range( 0, len(X) ):
    new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], new_T[i], new_C[i], new_CT[i] ) )
  new_file.close()



########################################3
#   OUTPUT 00
################################3

print int

fname= "dens_ob4.00.ct0.Ag.-1.n"
print fname
new_file= open( fname, "w" )
for i in range( 0, len(X) ):
  new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], n_T[i], n_C[i], n_CT[i] ) )
new_file.close()

fname= "dens_ob4.00.ct0.Ag.-1.p"
print fname
new_file= open( fname, "w" )
for i in range( 0, len(X) ):
  new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], p_T[i], p_C[i], p_CT[i] ) )
new_file.close()

fname= "dens_ob4.00.ct0.Ag.-1"
print fname
new_file= open( fname, "w" )
for i in range( 0, len(X) ):
  new_file.write("{:f}\t{:f}\t{:f}\t{:f}\n".format( X[i], p_T[i]+n_T[i], p_C[i]+n_C[i], p_CT[i]+n_CT[i] ) )
new_file.close()



