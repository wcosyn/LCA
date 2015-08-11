#!/usr/bin/python
# use totalint file width1 width2 width3 width4
import sys
import numpy as np

func = lambda x, s, max: max* np.sqrt(2/np.pi)*x*x / s**3 * np.exp( -(x)**2 / (2*s**2) ) # Target function

width1= float( sys.argv[3] )
width2= float( sys.argv[4] )
width3= float( sys.argv[5] )
width4= float( sys.argv[6] )
max = float( sys.argv[1] )

f = open( sys.argv[2] , 'w' )
f.write("# width 1: %(w1)f \t width 2: %(w2)f \t width 3: %(w3)f \t width 4: %(w4)f  [fm^-1]\n" %{"w1": width1, "w2":width2, "w3":width3, "w4":width4} )

for p in np.arange(0, 3, 0.01 ):
	f.write("%(p)f \t %(s1)f \t %(s2)f \t %(s3)f \t %(s4)f \n" %{"p": p, "s1": func(p, width1, max), "s2": func(p, width2, max), "s3": func(p, width3, max), "s4": func(p, width4, max) } )

f.close()

