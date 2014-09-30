#!/usr/bin/python
import sys
import numpy as np
import pylab as pl

file = open(sys.argv[1])

s1= [0]
p1= [0]
d1= [0]
f1= [0]
g1= [0]
h1= [0]
i1= [0]
matrix1= []
matrix1.append(s1)
matrix1.append(p1)
matrix1.append(d1)
matrix1.append(f1)
matrix1.append(g1)
matrix1.append(h1)
matrix1.append(i1)
s2= [0]
p2= [0]
d2= [0]
f2= [0]
g2= [0]
h2= [0]
i2= [0]
matrix2= []
matrix2.append(s2)
matrix2.append(p2)
matrix2.append(d2)
matrix2.append(f2)
matrix2.append(g2)
matrix2.append(h2)
matrix2.append(i2)
colLabels=( 'S', 'P', 'D', 'F', 'G', 'H', 'I+' )

rows1= 0
rows2= 0
number1= 0
number2= 0

norm1= np.array( [0.0]*7)
norm2= np.array( [0.0]*7)
while 1:
  lines = file.readlines(100000)
  if not lines:
    break
  for line in lines:
    if line[0] == '#':
      continue
    array= line.split('\t')

    l1= int( array[1] )
    number1= max( number1, l1+1)
    for col in range(9, len(array) ):
      val= float( array[col] )
      l = col-9
      if( l > 6 ):
        l= 6
      lvector= matrix1.pop( l )
      norm1[l1]+= val
      while len( lvector ) < l1+1:
	lvector.append(0)
      lvector[l1]+= val
      matrix1.insert( l, lvector )
      rows1= max( rows1, l+1 )

    l2= int( array[5] )
    number2= max( number2, l2+1)
    for col in range(9, len(array) ):
      val= float( array[col] )
      l = col-9
      if( l > 6 ):
        l= 6
      lvector= matrix2.pop( l )
      while len( lvector ) < l2+1:
	lvector.append(0)
      lvector[l2]+= val
      norm2[l2]+= val
      matrix2.insert( l, lvector )
      rows2= max( rows2, l+1 )

#print matrix

norms1 = np.array( [0.0]*number1 )
for row in range( 0, rows1 ):
  norms1+=matrix1[row]

norms2 = np.array( [0.0]*number2 )
for row in range( 0, rows2 ):
  norms2+=matrix2[row]

colours = [ 'r', 'b', 'g', 'y', 'c', 'm', 'k', 'w' ]
ind1 = np.arange( 0.75 , number1+0.75 , 1 )
ind2 = np.arange( 0.75 , number2+0.75 , 1 )
yoff1 = np.array( [0.0]*number1)
yoff2 = np.array( [0.0]*number2)
width = 0.5

fig1 = pl.figure(figsize=(18,6))
fig1.subplots_adjust(bottom=0.10, top=0.97, right=0.99, left=0.09, hspace=0.3)
ax1 = fig1.add_subplot(121)
ax1.set_axisbelow(True)
ax1.set_ylabel(' PN Pairs', fontsize=20)
ax1.set_xlabel( sys.argv[2] + ' shell', fontsize=20)

for row in range( 0, rows1 ):
  print row, matrix1[row]
  ax1.bar( ind1, matrix1[row]/norms1, width, bottom = yoff1, color = colours[row], label=colLabels[row])
  yoff1 = yoff1 + matrix1[row]/norms1
pl.ylim([0,1])
locs, labels= pl.xticks( [1,2,3,4,5,6,7], ['s', 'p', 'd', 'f', 'g', 'h', 'i'] )
pl.setp( labels, fontsize=20 )
pl.xlim( [-0.5, 0.5+max(number1, number2) ] )
pl.legend( loc= 'center left' )


ax2 = fig1.add_subplot(122)
ax2.set_axisbelow(True)
ax2.set_ylabel(' PN Pairs', fontsize=20)
ax2.set_xlabel( sys.argv[3] + ' shell', fontsize=20)

for row in range( 0, rows2 ):
  print row, matrix2[row]
  ax2.bar( ind2, matrix2[row]/norms2, width, bottom = yoff2, color = colours[row], label=colLabels[row])
  yoff2 = yoff2 + matrix2[row]/norms2

pl.ylim([0,1])
locs, labels= pl.xticks( [1,2,3,4,5,6,7], ['s', 'p', 'd', 'f', 'g', 'h', 'i'] )
pl.setp( labels, fontsize=20 )
pl.xlim( [-0.5, 0.5+max(number1, number2) ] )

pl.show()

