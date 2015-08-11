#!/usr/bin/env python
import sys
import numpy as np

def getCorrection( A ):

  result= 0
  #if( A <= 12 ):
  #  rico= (1.30-1)/10
  #  result= 1+ rico*(A-2)
  if( A <= 56 ):
    rico= (1.70-1.64)/44
    result= 1.64+ rico* (A-12)
#    rico= (1.39-1.30)/44
#    result= 1.30 + rico* (A-12)
  else:
    rico= (1.71-1.70)/152
    result= 1.70 + rico*(A-56)
  return result

def getCorrectionLow( A ):
  if ( A <= 56 ):
    rico= (0.27-0.23)/44
    result= 0.23+ rico* (A-12)
    return getCorrection( A ) - result
  else:
    rico= (0.29-0.27)/152
    result= 0.23+ rico* (A-56)
    return getCorrection( A ) - result

def getCorrectionHigh( A ):
  if ( A <= 56 ):
    rico= (0.27-0.23)/44
    result= 0.23+ rico* (A-12)
    return getCorrection( A ) + result
  else:
    rico= (0.29-0.27)/152
    result= 0.23+ rico* (A-56)
    return getCorrection( A ) + result




#file = sys.argv[1]
file = "PairsNP"
print file
data = np.loadtxt(file) 

A = data[:,0]
Z = data[:,1]
N = A - Z
# (n l S): 3=(- 0 -), 4=(- 0 0), 5=(- 0 1), 6=(0 0 1), 7=(n l S L)=(0 0 1 0) 8=(- 0 1 0)
#F = data[:,5]
F = data[:,sys.argv[1]]


result =  2.*F/A #* float(sys.argv[2])

#Now add CM correction factor

#resultlow= result*getCorrection(A)
#resulthigh= result*getCorrection(A)

outfile = open( file + '.a2', 'w')

outfile.write( "# A \t\t Z \t\t a_2((nlS)=(xxx)) \t\t a_2((nlS)(xxx))_low \t a_2((nlS)=(xxx))_high\n")
for i in range(0, len(result)):
	outfile.write( `A[i]` + ' \t\t' + `Z[i]` + ' \t\t' + `result[i]` + ' \t\t' + `result[i]*getCorrectionLow( A[i] )` + '\t\t' + `result[i]*getCorrectionHigh( A[i] )` + '\n')

outfile.close()
	
