#calculate Moshinsky Brackets
from __future__ import division
import math as ma
import numpy as np
from sympy import Integer, pi, sqrt, sympify
#from sage.rings.complex_number import ComplexNumber
#from sage.rings.finite_rings.integer_mod import Mod

# This list of precomputed factorials is needed to massively
# accelerate future calculations of the various coefficients
_Factlist = [1]


def _calc_factlist(nn):
    r"""
    Function calculates a list of precomputed factorials in order to
    massively accelerate future calculations of the various
    coefficients.

    INPUT:

    -  ``nn`` -  integer, highest factorial to be computed

    OUTPUT:

    list of integers -- the list of precomputed factorials

    EXAMPLES:

    Calculate list of factorials::

        sage: from sage.functions.wigner import _calc_factlist
        sage: _calc_factlist(10)
        [1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800]
    """
    if nn >= len(_Factlist):
        for ii in range(len(_Factlist), int(nn + 1)):
            _Factlist.append(_Factlist[ii - 1] * ii)
    return _Factlist[:int(nn) + 1]


def wigner_3j(j_1, j_2, j_3, m_1, m_2, m_3):
    r"""
    Calculate the Wigner 3j symbol `Wigner3j(j_1,j_2,j_3,m_1,m_2,m_3)`.

    INPUT:

    -  ``j_1``, ``j_2``, ``j_3``, ``m_1``, ``m_2``, ``m_3`` - integer or half integer

    OUTPUT:

    Rational number times the square root of a rational number.

    """
    
    if int(j_1 * 2) != j_1 * 2 or int(j_2 * 2) != j_2 * 2 or \
            int(j_3 * 2) != j_3 * 2:
        raise ValueError("j values must be integer or half integer")
    if int(m_1 * 2) != m_1 * 2 or int(m_2 * 2) != m_2 * 2 or \
            int(m_3 * 2) != m_3 * 2:
        raise ValueError("m values must be integer or half integer")
    if m_1 + m_2 + m_3 != 0:
        return 0
    prefid = Integer((-1) ** int(j_1 - j_2 - m_3))
    m_3 = -m_3
    a1 = j_1 + j_2 - j_3
    if a1 < 0:
        return 0
    a2 = j_1 - j_2 + j_3
    if a2 < 0:
        return 0
    a3 = -j_1 + j_2 + j_3
    if a3 < 0:
        return 0
    if (abs(m_1) > j_1) or (abs(m_2) > j_2) or (abs(m_3) > j_3):
        return 0

    maxfact = max(j_1 + j_2 + j_3 + 1, j_1 + abs(m_1), j_2 + abs(m_2),
                  j_3 + abs(m_3))
    _calc_factlist(int(maxfact))

    argsqrt = Integer(_Factlist[int(j_1 + j_2 - j_3)] *
                     _Factlist[int(j_1 - j_2 + j_3)] *
                     _Factlist[int(-j_1 + j_2 + j_3)] *
                     _Factlist[int(j_1 - m_1)] *
                     _Factlist[int(j_1 + m_1)] *
                     _Factlist[int(j_2 - m_2)] *
                     _Factlist[int(j_2 + m_2)] *
                     _Factlist[int(j_3 - m_3)] *
                     _Factlist[int(j_3 + m_3)]) / \
        _Factlist[int(j_1 + j_2 + j_3 + 1)]

    ressqrt = sqrt(argsqrt)
    if ressqrt.is_complex:
        ressqrt = ressqrt.as_real_imag()[0]

    imin = max(-j_3 + j_1 + m_2, -j_3 + j_2 - m_1, 0)
    imax = min(j_2 + m_2, j_1 - m_1, j_1 + j_2 - j_3)
    sumres = 0
    for ii in range(int(imin), int(imax) + 1):
        den = _Factlist[ii] * \
            _Factlist[int(ii + j_3 - j_1 - m_2)] * \
            _Factlist[int(j_2 + m_2 - ii)] * \
            _Factlist[int(j_1 - ii - m_1)] * \
            _Factlist[int(ii + j_3 - j_2 + m_1)] * \
            _Factlist[int(j_1 + j_2 - j_3 - ii)]
        sumres = sumres + Integer((-1) ** ii) / den

    res = ressqrt * sumres * prefid
    return res


def clebsch_gordan(j_1, j_2, j_3, m_1, m_2, m_3):
    r"""
    Calculates the Clebsch-Gordan coefficient
    `\langle j_1 m_1 \; j_2 m_2 | j_3 m_3 \rangle`.

    The reference for this function is [Edmonds74]_.

    INPUT:

    -  ``j_1``, ``j_2``, ``j_3``, ``m_1``, ``m_2``, ``m_3`` - integer or half integer

    OUTPUT:

    Rational number times the square root of a rational number.

    EXAMPLES::

        >>> from sympy import S
        >>> from sympy.physics.wigner import clebsch_gordan
        >>> clebsch_gordan(S(3)/2, S(1)/2, 2, S(3)/2, S(1)/2, 2)
        1
        >>> clebsch_gordan(S(3)/2, S(1)/2, 1, S(3)/2, -S(1)/2, 1)
        sqrt(3)/2
        >>> clebsch_gordan(S(3)/2, S(1)/2, 1, -S(1)/2, S(1)/2, 0)
        -sqrt(2)/2

    NOTES:

    The Clebsch-Gordan coefficient will be evaluated via its relation
    to Wigner 3j symbols:

    .. math::

        \langle j_1 m_1 \; j_2 m_2 | j_3 m_3 \rangle
        =(-1)^{j_1-j_2+m_3} \sqrt{2j_3+1} \;
        Wigner3j(j_1,j_2,j_3,m_1,m_2,-m_3)

    See also the documentation on Wigner 3j symbols which exhibit much
    higher symmetry relations than the Clebsch-Gordan coefficient.

    AUTHORS:

    - Jens Rasch (2009-03-24): initial version
    """
    res = (-1) ** sympify(j_1 - j_2 + m_3) * sqrt(2 * j_3 + 1) * \
        wigner_3j(j_1, j_2, j_3, m_1, m_2, -m_3)
    return res


def _big_delta_coeff(aa, bb, cc, prec=None):
    r"""
    Calculates the Delta coefficient of the 3 angular momenta for
    Racah symbols. Also checks that the differences are of integer
    value.

    INPUT:

    -  ``aa`` - first angular momentum, integer or half integer

    -  ``bb`` - second angular momentum, integer or half integer

    -  ``cc`` - third angular momentum, integer or half integer

    -  ``prec`` - precision of the ``sqrt()`` calculation

    OUTPUT:

    double - Value of the Delta coefficient

    EXAMPLES::

        sage: from sage.functions.wigner import _big_delta_coeff
        sage: _big_delta_coeff(1,1,1)
        1/2*sqrt(1/6)
    """

    if int(aa + bb - cc) != (aa + bb - cc):
        raise ValueError("j values must be integer or half integer and fulfill the triangle relation")
    if int(aa + cc - bb) != (aa + cc - bb):
        raise ValueError("j values must be integer or half integer and fulfill the triangle relation")
    if int(bb + cc - aa) != (bb + cc - aa):
        raise ValueError("j values must be integer or half integer and fulfill the triangle relation")
    if (aa + bb - cc) < 0:
        return 0
    if (aa + cc - bb) < 0:
        return 0
    if (bb + cc - aa) < 0:
        return 0

    maxfact = max(aa + bb - cc, aa + cc - bb, bb + cc - aa, aa + bb + cc + 1)
    _calc_factlist(maxfact)

    argsqrt = Integer(_Factlist[int(aa + bb - cc)] *
                     _Factlist[int(aa + cc - bb)] *
                     _Factlist[int(bb + cc - aa)]) / \
        Integer(_Factlist[int(aa + bb + cc + 1)])

    ressqrt = sqrt(argsqrt)
    if prec:
        ressqrt = ressqrt.evalf(prec).as_real_imag()[0]
    return ressqrt


def racah(aa, bb, cc, dd, ee, ff, prec=None):
    """
    Calculate the Racah symbol `W(a,b,c,d;e,f)`.

    INPUT:

    -  ``a``, ..., ``f`` - integer or half integer

    -  ``prec`` - precision, default: ``None``. Providing a precision can
       drastically speed up the calculation.

    OUTPUT:

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    Examples
    ========

    >>> from sympy.physics.wigner import racah
    >>> racah(3,3,3,3,3,3)
    -1/14

    NOTES:

    The Racah symbol is related to the Wigner 6j symbol:

    .. math::

       Wigner6j(j_1,j_2,j_3,j_4,j_5,j_6)
       =(-1)^{j_1+j_2+j_4+j_5} W(j_1,j_2,j_5,j_4,j_3,j_6)

    Please see the 6j symbol for its much richer symmetries and for
    additional properties.

    ALGORITHM:

    This function uses the algorithm of [Edmonds74]_ to calculate the
    value of the 6j symbol exactly. Note that the formula contains
    alternating sums over large factorials and is therefore unsuitable
    for finite precision arithmetic and only useful for a computer
    algebra system [Rasch03]_
    """
    prefac = _big_delta_coeff(aa, bb, ee, prec) * \
        _big_delta_coeff(cc, dd, ee, prec) * \
        _big_delta_coeff(aa, cc, ff, prec) * \
        _big_delta_coeff(bb, dd, ff, prec)
    if prefac == 0:
        return 0
    imin = max(aa + bb + ee, cc + dd + ee, aa + cc + ff, bb + dd + ff)
    imax = min(aa + bb + cc + dd, aa + dd + ee + ff, bb + cc + ee + ff)

    maxfact = max(imax + 1, aa + bb + cc + dd, aa + dd + ee + ff,
                 bb + cc + ee + ff)
    _calc_factlist(maxfact)

    sumres = 0
    for kk in range(int(imin), int(imax) + 1):
        den = _Factlist[int(kk - aa - bb - ee)] * \
            _Factlist[int(kk - cc - dd - ee)] * \
            _Factlist[int(kk - aa - cc - ff)] * \
            _Factlist[int(kk - bb - dd - ff)] * \
            _Factlist[int(aa + bb + cc + dd - kk)] * \
            _Factlist[int(aa + dd + ee + ff - kk)] * \
            _Factlist[int(bb + cc + ee + ff - kk)]
        sumres = sumres + Integer((-1) ** kk * _Factlist[kk + 1]) / den

    res = prefac * sumres * (-1) ** int(aa + bb + cc + dd)
    return res




def A(l1,l,l2,L,x):
     s = 0
     c = ((ma.factorial(l1+l+x+1)*ma.factorial(l1+l-x)*ma.factorial(l1+x-l))/ma.factorial(l+x-l1))**(0.5)
     d = ((ma.factorial(l2+L+x+1)*ma.factorial(l2+L-x)*ma.factorial(l2+x-L))/ma.factorial(L+x-l2))**(0.5)
     q = 0
     while((l+l1-q) >= 0 and (L+l2-q) >= 0):
          if((q-x) >= 0 and (L+l2-q) % 2 == 0 and (l+l1-q) % 2 == 0 ):
               s = s + (-1)**(0.5*(l+q-l1))*(ma.factorial(l+q-l1)/(ma.factorial((l+q-l1)/2)*ma.factorial((l+l1-q)/2)))*(1/(ma.factorial(q-x)*ma.factorial(q+x+1)))*(ma.factorial(L+q-l2)/(ma.factorial((L+q-l2)/2)*ma.factorial((L+l2-q)/2)))
          q = q + 1
     return c*d*s
     

     
def Coef(n, l, N, L, l1, l2):
     return ((ma.factorial(l1)*ma.factorial(l2)/(ma.factorial(2*l1)*ma.factorial(2*l2)))*((2*l+1)*(2*L+1)/2**(l+L))*(ma.factorial(n+l)/(ma.factorial(n)*ma.factorial(2*n+2*l+1)))*(ma.factorial(N+L)/(ma.factorial(N)*ma.factorial(2*N+2*L+1))))



def start(n, l, N, L, l1, l2, c):
     if(ma.fabs(l-l1)<= ma.fabs(L-l2)):
          x = ma.fabs(L-l2)
     else:
          x = ma.fabs(l-l1)
     if(l+l1 >= L+l2):
          y = L + l2
     else:
          y = l + l1 
     s = 0
     while(x <= y):
          s = s + ((2*x+1)*A(l1,l,l2,L,x)*racah(l,L,l1,l2,c,x))
          x = x + 1
     return (ma.sqrt(Coef(n,l,N,L,l1,l2)))*((-1)**(n+l+L-c))*s

def M_element(n, l, N, L, q, v, Q, V, c ):
     s = 0
     if(q == n-1 and v == l and Q == N and V == L):
          s = 0.5*ma.sqrt(n*(n+l+0.5))
     elif(q == n and v == l and Q == N-1 and V == L):
          s = 0.5*ma.sqrt(N*(N+L+0.5))
     elif(q == n-1 and v == l+1 and Q == N-1 and V == L+1):
          s = ma.sqrt(n*N*(l+1)*(L+1))*(-1)**(c+L+l)*racah(l,l+1,L,L+1,1,c)
     elif(q == n-1 and v == l+1 and Q == N and V == L-1):
          s = ma.sqrt(n*(N+L+0.5)*(l+1)*L)*(-1)**(c+L+l)*racah(l,l+1,L,L-1,1,c)
     elif(q == n and v == l-1 and Q == N-1 and V == L+1):
          s = ma.sqrt((n+l+0.5)*N*l*(L+1))*(-1)**(c+L+l)*racah(l,l-1,L,L+1,1,c)
     elif(q == n and v == l-1 and Q == N and V == L-1):
          s = ma.sqrt((n+l+0.5)*(N+L+0.5)*l*L)*(-1)**(c+L+l)*racah(l,l-1,L,L-1,1,c)
     else: 
          s = 0
     return s

def Moshinsky(n, l, N, L, n1, l1, n2, l2, c):
     if(2*n1 + l1 + 2*n2 + l2 == 2*n + l + 2*N + L):
          if(n1 == 0 and n2 == 0):
               return start(n, l, N, L, l1, l2, c)
          else:
               if(n1 > 0):
                    res = 0
                    if(n >= 1):
                         i = n-1
                    else:
                         i = n
                    while(i <= n):
                         if(l >= 1):
                              j = l-1
                         else:
                              j = l 
                         while(j <= l+1):
                              if(N >= 1):
                                   k = N-1
                              else:
                                   k = N
                              while(k <= N):
                                   if(L >= 1):
                                        m = L-1
                                   else:
                                        m = L
                                   while(m <= L+1):
                                        res = res + ((n1)*(n1+l1+1./2.))**(-0.5)*M_element(n,l,N,L,i,j,k,m,c)*Moshinsky(i,j,k,m,n1-1,l1,n2,l2,c)
                                        m = m + 1
                                   k = k + 1
                              j = j + 1
                         i = i + 1
               elif(n2 > 0 and n1 == 0):
                    res = 0
                    if(n >= 1):
                         i = n-1
                    else:
                         i = n
                    while(i <= n):
                         if(l >= 1):
                              j = l-1
                         else:
                              j = l 
                         while(j <= l+1):
                              if(N >= 1):
                                   k = N-1
                              else:
                                   k = N
                              while(k <= N):
                                   if(L >= 1):
                                        m = L-1
                                   else:
                                        m = L
                                   while(m <= L+1):
                                        res = res + ((n2)*(n2+l2+1./2.))**(-0.5)*(-1)**(L-c)*M_element(n,l,N,L,i,j,k,m,c)*Moshinsky(i,j,k,m,n2-1,l2,0,l1,c)
                                        m = m + 1
                                   k = k + 1
                              j = j + 1
                         i = i + 1
               return res

     else:
          return 0

print clebsch_gordan(2,1,3,1,1,2)
          
def format(value):
    return "%.9f" % value
             
n1 = np.array([0,0,0,0,0,0,2,2,2,2])
l1 = np.array([0,1,1,2,2,2,2,2,2,2])
n2 = np.array([0,0,0,0,0,0,1,1,1,1])
l2 = np.array([0,3,5,2,4,5,3,3,4,4])
N  = np.array([0,0,0,0,1,0,0,1,0,3])
L  = np.array([0,2,1,1,3,5,3,0,2,2])
n  = np.array([0,1,0,0,0,0,1,2,4,0])
l  = np.array([0,0,5,3,1,2,6,5,2,4])
lam= np.array([0,2,6,4,3,4,4,5,2,4])

f = open("moshinsky.txt", "w")
f.write("n" + '\t' + "l_1" + '\t' + "n_2" '\t' +  "l_2" + '\t' + "N" + '\t' + "L" +'\t' + "n" + '\t' + "l" +'\t' + "Lambda" + '\n')

for i in range(0,10):
     f.write(str(n1[i]) + '\t' + str(l1[i]) + '\t' + str(n2[i]) + '\t' +  str(l2[i]) + '\t' + str(N[i]) + '\t' + str(L[i]) +'\t' + str(n[i]) + '\t' + str(l[i]) +'\t' + str(lam[i]) + '\t' )
     f.write( format(Moshinsky(n[i],l[i],N[i],L[i],n1[i],l1[i],n2[i],l2[i],lam[i]))) 
     f.write('\n')

f.close()

