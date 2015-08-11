#bepalen van aantal nucleonen in elke nl-toestand
import math
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad

########################
## nuclide parameters ##
########################

nuc = ["He","Be","C","O","Al","Ar","Ca40","Ca48","Fe","Ag","Pb"]
A   = [4,9,12,16,27,40,40,48,56,108,208]
Z   = [2,4,6,8,13,18,20,20,26,47,82]
N = [m - n for m,n in zip(A,Z)] #number of neutrons

degen   = np.array([2,4,2,6,2,4,8,4,6,2,10,8,6,4,2,12,10,8,6,4,2,14,10,6,12,8,2,4])
orbital = np.array([0,1,1,2,0,2,3,1,3,1,4,4,2,2,0,5,5,3,3,1,1,6,4,2,6,4,0,2])
radial = np.array([0,0,0,0,1,0,0,1,0,1,0,0,1,1,2,0,0,1,1,2,2,0,1,2,0,1,3,2])

nos = []
for i in xrange(0,len(degen)):
     s = 0
     for j in xrange(0,i+1):
           s = s + degen[j]
     nos.append(s)
     
def level(Q):
     j = 0
     s = 0
     while(s < Q):
          s = s + degen[j]
          j = j + 1
     return j
     
def vulOp(i,t):
     opv =  [[0 for x in range(15)] for x in range(15)] 
     result = 0
     if(t == 0):
          Q = N[i]
     elif(t == 1):
          Q = Z[i]
     for lvl in xrange(0, level(Q)):
          n = radial[lvl]
          l = orbital[lvl]
          if(Q < nos[lvl]):
               opv[n][l] = opv[n][l] + Q-nos[lvl-1]
          else:
               opv[n][l] = opv[n][l] + degen[lvl]
     return opv