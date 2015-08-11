#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color "Helvetica, 29"
set output "hodistrPb.eps"

set key right top height 1
set border lw 2
set tics scale 1,1
set pointsize 2.5

set xlabel "P_{cm} (GeV)"  font "Times-Italic, 30"
set ylabel "P_{cm}^2 P_2(P_{cm}) (GeV^{-1})" font "Times-Italic, 30"
set yr[0:15000]

plot "hocommom.Pb" u 1:2 w l ti "^{208}Pb         " lw 4 , "hocommom.Pb" u 1:($3)*10 ti "^{208}Pb | ^1S_0 (x10)" lw 4 w l
