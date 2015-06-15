#!/usr/bin/python
# use totalint file kT1 kT2 kT3 
# kT in MeV
import sys
import numpy as np

Mn= 938. # MeV
hbarc= 197.327 # MeV

func = lambda x, kT, max: max* np.sqrt(2/np.pi)*x*x / Mn**1.5 / kT**1.5 * np.exp( -(x)**2 / (2*Mn*kT) ) # Target function

width1= float( sys.argv[3] )
width2= float( sys.argv[4] )
width3= float( sys.argv[5] )
max = float( sys.argv[1] )

f = open( sys.argv[2] , 'w' )
f.write("# width 1: %(w1)f \t width 2: %(w2)f \t width 3: %(w3)f  [fm^-1]\n" %{"w1": width1, "w2":width2, "w3":width3} )

# p in MeV
for p in np.arange(0*hbarc, 3*hbarc, 0.01*hbarc ):
	f.write("%(p)f \t %(s1)f \t %(s2)f \t %(s3)f \n" %{"p": p, "s1": func(p, width1, max), "s2": func(p, width2, max), "s3": func(p, width3, max) } )

f.close()

