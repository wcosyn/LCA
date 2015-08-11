#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color "Helvetica, 29"
set output "gausshoFe.eps"

set key right top height 1
set border lw 2
set tics scale 1,1
set pointsize 2.5

set xlabel "P_{cm} (GeV)"  font "Times-Italic, 30"
set ylabel "P_2(P_{cm} | ^1S_0) (GeV^{-3})" font "Times-Italic, 30"
set logscale y
set yr[1e1:1e5]


plot "hocommom.Fe" u  1:($3)/($1)/($1) ti "^{56}Fe | ^1S_0" lw 4 w l, \
"datahofit.Fe" u 1:($2)/($1)/($1):($3)/($1)/($1) with filledcurves lw 4  ti "{/Symbol s}= 180 +/- 4 "
