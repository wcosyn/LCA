#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color "Helvetica, 29"
set output "hodistrC.eps"

set key right top height 1 Left
set border lw 2
set tics scale 1,1
set pointsize 2.5

set xr[0:0.7]
set yr[0:200]
set xlabel "P_{12} (GeV)"  font "Times-Italic, 30"
set ylabel "P_{12}^2 P_2(P_{12}) (GeV^{-1})" font "Times-Italic, 30"

plot "hocommom.C" u 1:2 w l ti "^{12}C" lw 5 , "hocommom.C" u 1:($3)*5 ti "^{12}C | ^3S_1 (x5)" lw 5 w l
