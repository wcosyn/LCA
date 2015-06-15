#!/usr/bin/python
import sys
import numpy as np
import pylab as pl

file = open(sys.argv[1])

s= [0]
p= [0]
d= [0]
f= [0]
g= [0]
h= [0]
i= [0]



matrix= []
matrix.append(s)
matrix.append(p)
matrix.append(d)
matrix.append(f)
matrix.append(g)
matrix.append(h)
matrix.append(i)
colLabels=( 'S', 'P', 'D', 'F', 'G', 'H', 'I+' )
print matrix

rows= 0
number= 0
norm= np.array( [0.0]*7 )
while 1:
  lines = file.readlines(100000)
  if not lines:
    break
  for line in lines:
    if line[0] == '#':
      continue
    array= line.split('\t')
    l1= int( array[1] )
    l2= int( array[5] )
    number= max( number, l1+1, l2+1)
    for col in range(9, len(array) ):
      val= float( array[col] )
      l = col-9
      if( l > 6 ):
	l= 6
      lvector= matrix.pop( l )
      while len( lvector ) < l1+1:
	lvector.append(0)
      lvector[l1]+= val
      norm[l1]+= val
      if( l1 != l2 ):
	while len( lvector ) < l2+1:
	  lvector.append(0)
	lvector[l2]+= val
        norm[l2]+= val
      matrix.insert( l, lvector )
      rows= max( rows, l+1 )

print rows, number
#print matrix

colours = [ 'r', 'b', 'g', 'y', 'c', 'm', 'k', 'w' ]
ind = np.arange( 0.75 , number+0.75 , 1 )
yoff = np.array( [0.0]*number)
width = 0.5

fig = pl.figure(figsize=(9,6))
fig.subplots_adjust(bottom=0.14, top=0.97, right=0.99, left=0.09)
#fig.clf()
ax1 = fig.add_subplot(111)
ax1.set_axisbelow(True)
#ax1.set_title(title)
ax1.set_ylabel('Number of ' + sys.argv[2] + ' ' + sys.argv[2] + ' Pairs', fontsize=18)
ax1.set_xlabel( sys.argv[2] + ' shell' )

norms = np.array( [0.0]*number )
for row in range( 0, rows ):
  norms+=matrix[row]

print "norms", norms
print "norm", norm

for row in range( 0, rows ):
  print row, matrix[row]
  pl.bar( ind, matrix[row], width, bottom = yoff, color = colours[row], label=colLabels[row])
  yoff = yoff + matrix[row]

#pl.ylim( [0,1 ] )
locs, labels= pl.xticks( [1,2,3,4,5,6,7], ['s', 'p', 'd', 'f', 'g', 'h', 'i'] )
pl.setp( labels, fontsize=20 )
pl.xlim( [-0.5, 0.5+number ] )

pl.legend( loc= 'center left')
pl.show()

