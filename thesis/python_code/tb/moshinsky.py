#calculate Moshinsky Brackets
from __future__ import division
import math as ma
import numpy as np
from sympy import Integer, pi, sqrt, sympify
from sympy.physics import wigner
#from sage.rings.complex_number import ComplexNumber
#from sage.rings.finite_rings.integer_mod import Mod

# This list of precomputed factorials is needed to massively
# accelerate future calculations of the various coefficients
_Factlist = [1]


def _calc_factlist(nn):
    if nn >= len(_Factlist):
        for ii in range(len(_Factlist), int(nn + 1)):
            _Factlist.append(_Factlist[ii - 1] * ii)
    return _Factlist[:int(nn) + 1]


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
          s = s + ((2*x+1)*A(l1,l,l2,L,x)*wigner.racah(l,L,l1,l2,c,x))
          x = x + 1
     return (ma.sqrt(Coef(n,l,N,L,l1,l2)))*((-1)**(n+l+L-c))*s

def M_element(n, l, N, L, q, v, Q, V, c ):
     s = 0
     if(q == n-1 and v == l and Q == N and V == L):
          s = 0.5*ma.sqrt(n*(n+l+0.5))
     elif(q == n and v == l and Q == N-1 and V == L):
          s = 0.5*ma.sqrt(N*(N+L+0.5))
     elif(q == n-1 and v == l+1 and Q == N-1 and V == L+1):
          s = ma.sqrt(n*N*(l+1)*(L+1))*(-1)**(c+L+l)*wigner.racah(l,l+1,L,L+1,1,c)
     elif(q == n-1 and v == l+1 and Q == N and V == L-1):
          s = ma.sqrt(n*(N+L+0.5)*(l+1)*L)*(-1)**(c+L+l)*wigner.racah(l,l+1,L,L-1,1,c)
     elif(q == n and v == l-1 and Q == N-1 and V == L+1):
          s = ma.sqrt((n+l+0.5)*N*l*(L+1))*(-1)**(c+L+l)*wigner.racah(l,l-1,L,L+1,1,c)
     elif(q == n and v == l-1 and Q == N and V == L-1):
          s = ma.sqrt((n+l+0.5)*(N+L+0.5)*l*L)*(-1)**(c+L+l)*wigner.racah(l,l-1,L,L-1,1,c)
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
                                        if(i == n-1 and j == l and k == N and m == L):
                                             res = res + ((n2)*(n2+l2+1./2.))**(-0.5)*M_element(n,l,N,L,i,j,k,m,c)*Moshinsky(i,j,k,m,0,l1,n2-1,l2,c)
                                        elif(i == n and j == l and k == N-1 and m == L):
                                             res = res + ((n2)*(n2+l2+1./2.))**(-0.5)*M_element(n,l,N,L,i,j,k,m,c)*Moshinsky(i,j,k,m,0,l1,n2-1,l2,c)
                                        else:
                                             res = res + ((n2)*(n2+l2+1./2.))**(-0.5)*((-1)*M_element(n,l,N,L,i,j,k,m,c))*Moshinsky(i,j,k,m,0,l1,n2-1,l2,c)
                                        m = m + 1
                                   k = k + 1
                              j = j + 1
                         i = i + 1
               return res

     else:
          return 0


          
def format_(value):
    return "%.9f" % value

             
"""n1 = np.array([0,0,0,0,0,0,2,2,2,2])
l1 = np.array([0,1,1,2,2,2,2,2,2,2])
n2 = np.array([0,0,0,0,0,0,1,1,1,1])
l2 = np.array([0,3,5,2,4,5,3,3,4,4])
N  = np.array([0,0,0,0,1,0,0,1,0,3])
L  = np.array([0,2,1,1,3,5,3,0,2,2])
n  = np.array([0,1,0,0,0,0,1,2,4,0])
l  = np.array([0,0,5,3,1,2,6,5,2,4])
lam= np.array([0,2,6,4,3,4,4,5,2,4]) """


n1 = 1
l1 = 2
n2 = 2
l2 = 3
f = open("./mosh_test/recmosh" + str(n1)+ str(l1) + str(n2)+ str(l2) + ".txt", "w")
for n in xrange(0, 2*n1 + l1 + 2*n2 + l2 +1):
     for l in xrange(0, 2*n1 + l1 + 2*n2 + l2 +1 - 2*n):
          for N in xrange(0, 2*n1 + l1 + 2*n2 + l2 +1 -l - 2*n):
               for L in xrange(0, 2*n1 + l1 + 2*n2 + l2 +1 - l- 2*n-2*N):
                    print L
                    if(2*n1 + l1 + 2*n2 + l2  == 2*n+ l+ 2*N + L):
                         for lambd in xrange(abs(l1-l2), l1 + l2 +1):
                              f.write(str(100000*n +10000*l+1000*N+ 100*L + 10*lambd))
                              f.write('\t')
                              f.write(format_(Moshinsky(n,l,N,L,n1,l1,n2,l2,lambd)))
                              f.write('\n')

f.close()



