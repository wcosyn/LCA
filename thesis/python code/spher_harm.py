#spherical harmonics
import math as ma
import scipy as sp
import numpy as np
from scipy import special
from scipy.integrate import dblquad

l_1 = int(input("l: "))
m_1 = int(input("m: "))
l_2 = int(input("l': "))
m_2 = int(input("m': "))


def spherical_harmonic(theta, phi, l,m):
     return sp.special.sph_harm(m,l, phi, theta)

def integrand(theta, phi, l_1,m_1, l_2, m_2):
     return np.sin(theta)*spherical_harmonic(theta, phi, l_1,m_1)*spherical_harmonic(theta, phi, l_2,m_2)

norm, err  = dblquad(integrand, 0, 2*np.pi,
                   lambda x: 0,
                   lambda x: np.pi, args=( l_1,m_1, l_2, m_2))

print norm