import numpy as np

nuclide = int(input("Nuclide? (He, C, Al, Fe, Pb): "))

nuc = ["He","C","Al","Fe","Pb"]

def distr(nucl):
     result = 0
     for i in xrange(0,2):
          d = np.loadtxt("hocommom_" + "np" + "_" + "-1" + str(i) + "." + nuc[nucl])
          result =  result + d[:,1]
     return result
     
def NN(nucl):
     a = np.loadtxt("hocommom_" + "nn" + "_" + "-1-1"  + "." + nuc[nucl])
     b = np.loadtxt("hocommom_" + "pp" + "_" + "-1-1"  + "." + nuc[nucl]) 
     return a[:,1] + b[:,1]
               

d_array = distr(nuclide)  + NN(nuclide)
data = np.loadtxt("hocommom_" + "pp" + "_" + "-1-1"  + "." + nuc[nuclide]) 
k_array = data[:,0]

fn = nuc[nuclide] + "_cm" + "_maarten.txt"
f = open(fn,"w")
f.write("#k"+ '\t' + "cm_distr" + '\n')
for i in xrange(0,len(d_array)):
     f.write(str(k_array[i]) + '\t' +  str(d_array[i]) + '\n')
f.close()

     
     
