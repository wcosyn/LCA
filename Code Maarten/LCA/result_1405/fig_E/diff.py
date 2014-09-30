#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
import glob
import os


if len(sys.argv) == 2:

  A=sys.argv[1]

  all= np.loadtxt("./dens_rel.00.110.{:s}.-1-1-1-1.-1-1".format(A,A))
  X= all[:,0]
  allmf= all[:,1]
  allco2b= all[:,2]
  allco3b= all[:,3]
  allto= all[:,4]
  p00= np.loadtxt("./dens_rel.00.110.{:s}.0000.-1-1".format(A,A))
  mf00= p00[:,1]
  co2b00= p00[:,2]
  co3b00= p00[:,3]
  to00= p00[:,4]

  mfdiff= allmf-mf00
  co2bdiff= allco2b-co2b00
  co3bdiff= allco3b-co3b00
  todiff= allto-to00

  #newf= open("{:s}/dens_ob3.00.110.{:s}.-1".format(A,A), "w")
  newf= open("./dens_rel.00.110.{:s}.diff".format(A), "w")


  for i in range(0, len(X)):
    newf.write("%f\t%f\t%f\t%f\t%f\n" % ( X[i], mfdiff[i], co2bdiff[i], co3bdiff[i], todiff[i] ) )

  dX= X[1]-X[0]
  iall=dX*sum(X*X*(allto))
  i00=dX*sum(X*X*(to00))
  idiff=dX*sum(X*X*(todiff))
  print "integral: ", iall, i00, iall-i00
  print "diff integral: " , idiff

if len(sys.argv) >= 3:
  A=sys.argv[1]

  all= np.loadtxt("./dens_rel.00.110.{:s}.-1-1-1-1.-1-1".format(A,A))
  X= all[:,0]
  dX= X[1]-X[0]
  allmf= all[:,1]
  allco2b= all[:,2]
  allco3b= all[:,3]
  allto= all[:,4]

  p00= np.loadtxt("./dens_rel.00.110.{:s}.0000.-1-1".format(A,A))
  mf00= p00[:,1]
  co2b00= p00[:,2]
  co3b00= p00[:,3]
  to00= p00[:,4]

  mfdiff= allmf-mf00
  co2bdiff= allco2b-co2b00
  co3bdiff= allco3b-co3b00
  todiff= allto-to00

  i00=dX*sum(X*X*(co2b00+co3b00))
  for i in range( 2, len(sys.argv) ):
    print sys.argv[i]
    p= np.loadtxt(sys.argv[i])
    mf= p[:,1]
    co2b= p[:,2]
    co3b= p[:,3]
    to= p[:,4]
    mfdiff -= mf
    co2bdiff -= co2b
    co3bdiff -= co3b
    todiff -= to
    i00+=dX*sum(X*X*(co2b+ co3b))


  #newf= open("{:s}/dens_ob3.00.110.{:s}.-1".format(A,A), "w")
  newf= open("./dens_rel.00.110.{:s}.diff".format(A), "w")


  for i in range(0, len(X)):
    newf.write("%f\t%f\t%f\t%f\t%f\n" % ( X[i], mfdiff[i], co2bdiff[i], co3bdiff[i], todiff[i] ) )

  iall=dX*sum(X*X*(allco2b+allco3b))
  idiff=dX*sum(X*X*(co2bdiff+ co3bdiff))
  print "integral: ", iall, i00, iall-i00
  print "diff integral: " , idiff
