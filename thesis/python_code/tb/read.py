#read data files that contain moshinsky brackets (recmoshXXXX.dat)

import numpy as np

def choose_file( n_1, l_1, n_2, l_2):
     file_number = str(n_1) + str(l_1) + str(n_2) + str(l_2)
     file_name = "./mosh_b/recmosh" + file_number + ".dat"
     return np.loadtxt(file_name,skiprows=1)
#returns 2D array [row][100000*n + 10000*l + 1000*N + 100*l + 10*Lambda, Mosh_bracket]
#first element [0][1000*n_1 + 100*l_1 + 10*n_2 + 1*l_2]


def moshinsky( n_1, l_1, n_2, l_2, n, l, N, L, Lambd):
     result = 0
     if(n_1*10 + l_1 > n_2*10 + l_2):
          z = choose_file(n_2, l_2, n_1, l_1)
     else:
          z = choose_file(n_1, l_1, n_2, l_2)
     number = 100000*n + 10000*l + 1000*N + 100*L + 10*Lambd
     for i in xrange(0, len(z)):
          if(number == z[i][0]):
               if(abs(z[i][1]) > 0.000000000001):
                    if(n_1*10 + l_1 > n_2*10 + l_2):
                         result = ((-1.0)**(L-Lambd))*z[i][1]
                    else:
                         result = z[i][1]
               else:
                    result = 0
               break
     return result
 
f = open("mosh.txt","w")
f.write("l_1" + '\t' + "l_2" + '\t' + "lam" + '\t' + "n" + '\t' + "l"  + '\t'+ "N" + '\t' + "L")
f.write('\n')    
for l_1 in xrange(0,2):
     for l_2 in xrange(0,2):
          for n in xrange(0, l_1 + l_2 +1):
               for l in xrange(0, l_1 + l_2 +1):
                    for N in xrange(0, l_1 + l_2 +1):
                         for L in xrange(0, l_1 + l_2 +1):
                              if(l_1 + l_2 == 2*n+ l+ 2*N + L):
                                   for lambd in xrange(abs(l_1-l_2), l_1 + l_2 +1):
                                        f.write(str(l_1) + '\t'+ str(l_2) + '\t'+ str(lambd) + '\t'+ str(n) + '\t'+ str(l) + '\t'+ str(N) + '\t'+ str(L) + '\t'+ str(moshinsky(0, l_1, 0, l_2, n, l, N, L, lambd)) )
                                        f.write('\n')
f.close()
                              
               
               
               

"""f = open("mosh.txt","w")
f.write("l_1" + '\t' + "l_2" + '\t' + "lam" + '\t' + "n" + '\t' + "l"  + '\t'+ "N" + '\t' + "L")
f.write('\n')
n_1 = 2
n_2 = 1
for l_1 in xrange(0,5):
     for l_2 in xrange(0,5):
          for lambd in xrange(abs(l_1-l_2),l_1 + l_2 + 1):
               for n in xrange(0,2*n_1 + l_1 + 2*n_2 + l_2 +1 ):
                    for l in xrange(0,2*n_1 + l_1 + 2*n_2 + l_2 +1 -n  ):
                          for L in xrange(0,2*n_1 + l_1 + 2*n_2 + l_2 +1 -n  - l):
                               for N in xrange(0,2*n_1 + l_1 + 2*n_2 + l_2 + 1 -2*n-l-L ):
                                   if(2*n_1 + l_1 + 2*n_2 + l_2 == 2*n + l + 2*N + L):
                                        if(moshinsky(n_1, l_1, n_2, l_2, n, l, N, L, lambd) != 0):
                                             f.write(str(l_1) + '\t'+ str(l_2) + '\t'+ str(lambd) + '\t'+ str(n) + '\t'+ str(l) + '\t'+ str(N) + '\t'+ str(L) + '\t'+ str(moshinsky(n_1, l_1, n_2, l_2, n, l, N, L, lambd)) )
                                             f.write('\n')
f.close()"""
                                        
                                        
                                        
