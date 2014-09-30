#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color "Helvetica, 29"
set output "gausshoC.eps"

set key right top height 1
set border lw 2
set tics scale 1,1
set pointsize 2.5

set xlabel "P_{cm} (GeV)"  font "Times-Italic, 30"
set ylabel "P_2(P_{cm} | ^1S_0) (GeV^{-3})" font "Times-Italic, 30"
set logscale y
set yr[1e0:1e4]



plot "hocommom.C" u 1:($3)/($1)/($1) ti "^{12}C | ^1S_0" lw 4 w l, \
"datahofit.C" u 1:($2)/($1)/($1):($3)/($1)/($1) with filledcurves lw 4  ti "{/Symbol s}= 161 +/- 4 "
