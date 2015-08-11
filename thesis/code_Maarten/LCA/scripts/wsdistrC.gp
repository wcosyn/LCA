#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color "Helvetica, 29"
set output "wsdistrC.eps"

set key right top height 1
set border lw 2
set tics scale 1,1
set pointsize 2.5

set xlabel "P_{cm} (GeV)"  font "Times-Italic, 30"
set ylabel "P_{cm}^2 P_2(P_{cm}) (GeV^{-1})" font "Times-Italic, 30"

plot "wscommom.C" u 1:2 w l ti "^{12}C WS" lw 4 , "wscommom.C" u 1:3 ti "^{12}C | ^1S_0 WS" lw 4 w l
