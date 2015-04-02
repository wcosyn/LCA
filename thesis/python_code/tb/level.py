import math
import scipy as sp
import numpy as np
from scipy import special
from sympy.physics.quantum.cg import CG
from scipy.integrate import quad


nuc = ["He","Be","C","O","Al","Ca40","Ca48","Fe","Ag","Ar","Pb"]
A   = [4,9,12,16,27,40,48,56,109,40,208]
Z   = [2,4,6,8,13,20,20,26,47,18,82]
N = [m - n for m,n in zip(A,Z)] #number of neutrons

Omega = []
for i in range(len(A)):
     Omega.append((45.0*(A[i]**(-1./3.)))-(25.0*(A[i]**(-2./3.))))
     
hbarc = 197.327 #(MeV*fm)
M_n = 939.565 #mass neutron (MeV) [Particle Data Group]
M_p = 938.272 #mass proton (MeV) [Particle Data Group]
M_average = (M_n + M_p)/2.

degen   = [2,4,2,6,2,4,8,4,6,2,10,8,6,4,2,12,10,8,6,4,2,14,10,6,12,8,2,4]
orbital = [0,1,1,2,0,2,3,1,3,1,4,4,2,2,0,5,5,3,3,1,1,6,4,2,6,4,0,2]
radial = [0,0,0,0,1,0,0,1,0,1,0,0,1,1,2,0,0,1,1,2,2,0,1,2,0,1,3,2] 
nos = []
for i in xrange(0,len(degen)):
     s = 0
     for j in xrange(0,i+1):
           s = s + degen[j]
     nos.append(s)
print nos     
def level(Q):
     j = 0
     s = 0
     while(s < Q):
          s = s + degen[j]
          j = j + 1
     return j
       
"""k = open("tabel.txt", "w")  
k.write("#t_1" + '\t' +  "t_2" + "n_1" + '\t' + "l_1" +  '\t'+ "n_2" + '\t' + "l_2" + '\t' +"d_1" + '\t' +"d_2"+ '\t' +  "p? "+ '\t' +  "numbr of states" )
k.write('\n') """
def opvulling(i):
     result = 0
     for t_1 in xrange(0,2):
          if(t_1 == 0):
               Q_1 = N[i]
          elif(t_1 == 1):
               Q_1 = Z[i]
          for lvl_1 in xrange(0, level(Q_1)):
               n_1 = radial[lvl_1]
               l_1 = orbital[lvl_1]
               a = float(degen[lvl_1])/ float(2*(2*l_1 + 1))
               if(nos[lvl_1]-Q_1 > 0):
                    b = float(Q_1-nos[lvl_1-1])/float(degen[lvl_1])
               else:
                    b = 1
               for t_2 in xrange(0,2 ):
                    if(t_2 == 0):
                         Q_2 = N[i]
                    elif(t_2 == 1):
                         Q_2 = Z[i]
                    for lvl_2 in xrange(0,level(Q_2)):
                         n_2 = radial[lvl_2]
                         l_2 = orbital[lvl_2]
                         if(t_1 == t_2 and n_1 == n_2 and l_1 == l_2 and degen[lvl_2] == degen[lvl_1]):
                              c = float(degen[lvl_2]-1)/ ((2*(2*l_2 + 1))-1)
                         else:
                              c = float(degen[lvl_2])/(2*(2*l_2 + 1))
                         if(nos[lvl_2]-Q_2 > 0):
                              if(t_1 == t_2 and n_1 == n_2 and l_1 == l_2 and degen[lvl_2] == degen[lvl_1]):
                                   d = float(Q_2-nos[lvl_2-1]-1)/float(degen[lvl_2]-1)
                              else:
                                   d = float(Q_2-nos[lvl_2-1])/float(degen[lvl_2])
                         else:
                              d = 1   
                         for m_1 in xrange(-l_1, l_1 + 1):
                              for m_2 in xrange(-l_2, l_2 + 1):
                                   for s_1 in xrange(0,2):
                                        for s_2 in xrange(0,2):
                                             if(t_1 != t_2 or n_1 != n_2 or l_1 != l_2 or m_1 != m_2 or s_1 != s_2 or degen[lvl_2] != degen[lvl_1]):
                                                  """if(t_1 != t_2 or n_1 != n_2 or l_1 != l_2):
                                                  p = "no"
                                                  else:
                                                  p = "yes" """
                                                  """if(t_1 != t_2 or n_1 != n_2 or l_1 != l_2 or s_1 != s_2 or m_1 != m_2):"""
                                                  result = result + a*b*c*d 
                                               
                         """txt = str(t_1) +'\t' + str(t_2) +'\t'+ str(n_1) + '\t' + str(l_1)  + '\t'+ str(n_2) + '\t' + str(l_2)  + '\t'+ str(degen[lvl_1]) + '\t'+ str(degen[lvl_2]) + '\t' +str(result)+ '\t' + str(c*d*a*b)
                         k.write(txt)
                         k.write('\n')
                         k.write('\n')"""
     return result


#k.close()

for i in range(0,len(A)):
     print str(opvulling(i)) + '\t' + str(A[i]*(A[i]-1)) +'\t' +nuc[i]
"""som = 0    
z = np.loadtxt("tabel.txt")
for i in xrange(0,len(z)):
     som = som + z[i][8]
print som """
     
   
