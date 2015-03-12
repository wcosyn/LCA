import math
import scipy as sp
import numpy as np
from sympy.physics.quantum.cg import CG
from sympy import S



# opvulling van matrix met glebsch gordan coefficienten voor 2 spin-1/2 deeltjes clebsch_gordan_spin[S][ms_1+ 0.5][ms_2 + 0.5]
clebsch_gordan_spin = np.zeros((2,2,2))
for i in xrange(0,2):
	for j in xrange(0,2):
		for k in xrange(0,2):
			cg = CG(S(1)/2, S(j-0.5), S(1)/2, S(k-0.5), S(i), S(j+k-1))
			clebsch_gordan_spin[i][j][k] = cg.doit()

# opvulling van matrix met glebsch gordan coefficienten voor 2 deeltjes met arbitraire gehele l (conventie: clebsch_gordan[L][l_1][l_2][m_1 + 6][m_2 + 6] , M_L = m_1 + m_2 ligt dan vast)			
clebsch_gordan = np.zeros((13,13,13,13,13))
for i in xrange(0,13):
	for j in xrange(0,7):
		for k in xrange(0,7):
			for l in xrange(0,13):
				for m in xrange(0,13):
					cg = CG(j, l-6, k, m-6, i, m + l - 12)
					clebsch_gordan[i][j][k][l][m] = cg.doit()




