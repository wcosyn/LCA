#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import glob
import os


if len(sys.argv) == 2:

  A=sys.argv[1]

  all= np.loadtxt("./dens_ob4.00.111.{:s}.-1-1-1-1".format(A,A))
  X= all[:,0]
  allmf= all[:,1]
  allco= all[:,2]
  allto= all[:,3]
  p00= np.loadtxt("./dens_ob4.00.111.{:s}.0000".format(A,A))
  mf00= p00[:,1]
  co00= p00[:,2]
  to00= p00[:,3]

  mfdiff= allmf-mf00
  codiff= allco-co00
  todiff= allto-to00

  #newf= open("{:s}/dens_ob3.00.110.{:s}.-1".format(A,A), "w")
  newf= open("./dens_ob4.00.111.{:s}.diff".format(A), "w")


  for i in range(0, len(X)):
    newf.write("%f\t%f\t%f\t%f\n" % ( X[i], mfdiff[i], codiff[i], todiff[i] ) )

  dX= X[1]-X[0]
  iall=dX*sum(X*X*(allto))
  i00=dX*sum(X*X*(to00))
  idiff=dX*sum(X*X*(todiff))
  print "integral: ", iall, i00, iall-i00
  print "diff integral: " , idiff

if len(sys.argv) >= 3:
  A=sys.argv[1]

  all= np.loadtxt("./dens_ob4.00.111.{:s}.-1-1-1-1".format(A,A))
  X= all[:,0]
  dX= X[1]-X[0]
  allmf= all[:,1]
  allco= all[:,2]
  allto= all[:,3]

  p00= np.loadtxt("./dens_ob4.00.111.{:s}.0000".format(A,A))
  mf00= p00[:,1]
  co00= p00[:,2]
  to00= p00[:,3]

  mfdiff= allmf-mf00
  codiff= allco-co00
  todiff= allto-to00

  i00=dX*sum(X*X*(to00))
  for i in range( 2, len(sys.argv) ):
    print sys.argv[i]
    p= np.loadtxt(sys.argv[i])
    mf= p[:,1]
    co= p[:,2]
    to= p[:,3]
    mfdiff -= mf
    codiff -= co
    todiff -= to
    i00+=dX*sum(X*X*(to))


  #newf= open("{:s}/dens_ob3.00.110.{:s}.-1".format(A,A), "w")
  newf= open("./dens_ob4.00.111.{:s}.diff".format(A), "w")


  for i in range(0, len(X)):
    newf.write("%f\t%f\t%f\t%f\n" % ( X[i], mfdiff[i], codiff[i], todiff[i] ) )

  iall=dX*sum(X*X*(allto))
  idiff=dX*sum(X*X*(todiff))
  print "integral: ", iall, i00, iall-i00
  print "diff integral: " , idiff
