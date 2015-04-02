import math
import scipy as sp
import numpy as np
from scipy import special
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_6j, wigner_9j

     

def six_j(j_1, j_2, j_3, j_4, j_5, j_6):
     result = 0
     if( j_1 >= abs(j_2 - j_3) and j_1 <= j_2 + j_3 and j_1 >= abs(j_5 - j_6) and j_1 <= j_5 + j_6 and j_4 >= abs(j_2 - j_6) and j_4 <= j_2 + j_6 and j_4 >= abs(j_5 - j_3) and j_4 <= j_5 + j_3 ):
          if((j_1 + j_2 + j_3) % 1 == 0 and (j_1 + j_5 + j_6) % 1 == 0 and (j_4 + j_2 + j_6) % 1 == 0 and (j_4 + j_5 + j_3) % 1 == 0 ):
               result = wigner_6j(j_1, j_2, j_3, j_4, j_5, j_6,prec=10)
     return result
     
def nine_j(j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9):
     result = 0
     if(abs(j_1 - j_2) <= j_3 and j_1 + j_2 >= j_3 and abs(j_4 - j_5) <= j_6 and j_4 + j_5 >= j_6 and abs(j_7 - j_8) <= j_9 and j_7 + j_8 >= j_9 ):
          if(abs(j_1 - j_4) <= j_7 and j_1 + j_4 >= j_7 and abs(j_2 - j_5) <= j_8 and j_2 + j_5 >= j_8 and abs(j_3 - j_6) <= j_9 and j_3 + j_6 >= j_9 ):
               if(((j_1 + j_2 + j_3) % 1) == 0 and ((j_4 + j_5 + j_6) % 1) == 0 and ((j_7 + j_8 + j_9) % 1) == 0):
                    result = wigner_9j(j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9,prec=10)
               else: 
                    result = 0
          else: 
               result = 0
     else:
          result = 0
     return result
     
print nine_j(0.5, 2.5, 2, 1.5, 1, 1.5, 2, 1.5, 2.5)
     
