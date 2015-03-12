#read data files that contain moshinsky brackets (recmoshXXXX.dat)

import numpy as np

def choose_file(n_1,l_1,n_2,l_2):
     file_number = str(n_1) + str(l_1) + str(n_2) + str(l_2)
     file_name = "recmosh" + file_number + ".dat"
     return np.loadtxt(file_name)
     
#returns 2D array [row][n,l,N,l,Lambda,0, Mosh_bracket]
for l_1 in xrange(0,2):
     for l_2 in xrange(0,2):
          z = choose_file(0,l_1,0,l_2)
          for n in xrange(0,2):
               for l in xrange(0,3):
                    for N in xrange(0,2):
                         for L in xrange(0,3):
                              for lambd in xrange(0,5):
                                   for i in xrange(0,len(z)):
                                        if(n == z[i][0] and l == z[i][1] and N == z[i][2] and L == z[i][3] and lambd == z[i][4]):
                                             if(abs(z[i][6]) > 0.0000001):
                                                  mosh_bracket = z[i][6]
                                                  print "l_1 = " + str(l_1)
                                                  print "l_2 = " + str(l_2)
                                                  print "n = " + str(n)
                                                  print "l = " + str(l)
                                                  print "N = " + str(N)
                                                  print "L = " + str(L)
                                                  print "lambda = " + str(lambd)
                                                  print mosh_bracket
                                                  print "############"
