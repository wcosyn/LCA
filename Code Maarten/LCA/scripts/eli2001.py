#!/usr/bin/python

import sys
import numpy as np

data = np.loadtxt( sys.argv[1] )


f = open( "test.dat" , 'w' )
f.write( "# Centre Of Mass Momentum Distribution Of %(A)s Proton-Proton Pairs In R = sqrt(X**2 + Y**2 ) Direction \n"% {"A": sys.argv[2] } )

x = data[:, 0]
y = data[:, 2]
z = data[:, 1]

integral = 0
integralt = 0
for i in range( len(x) ):
	if i == 0:
		dx = x[1]-x[0]
		integral+= x[1]*dx
		integralt += z[1]* dx
	else :
		dx = x[i] - x[i-1]
		integral += y[i]* dx
		integralt += z[i]* dx
#f.write( "# 1S0 Distribution is Normalized to  %(i)f \n# = Number of 1S0 - Pairs.\n# Total Distribution is Normalized to Total Number of Proton-Proton Pairs: Z(Z-1)/2 \n"% { "i": integral } )
f.write( "# Model: %(model)s \n" %{ "model" : sys.argv[3] } )
f.write( "# C.O.M. Momentum ( GeV ) \t 1S0 Distribution ( GeV ** -2 )\t Total Distribution( GeV ** -2)\n" )

print "integral is ", integral, integralt


for i in range( len(x) ):
	if( i == 0 ):
		result = 0
		rest= 0
	else :
		result = y[i]/x[i]
		rest=  z[i]/x[i]
	#result = y[i]/x[i]/x[i]
	#rest=  z[i]/x[i]/x[i]
	#f.write("%(x)f \t\t\t %(y)f \n" %{ "x": x[i], "y": result} )
	f.write("%(x)f \t\t\t %(y)f \t\t\t %(z)f \n" %{ "x": x[i], "y": result, "z" : rest } )

f.close()

