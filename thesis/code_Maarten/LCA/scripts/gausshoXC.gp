#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color "Helvetica, 29"
set output "gaussXhoC.eps"

set style data line 
set key right top height 1
set border lw 2
set tics scale 1,1
set pointsize 2.5
set key spacing 2


set xlabel "P_{cm} (GeV)"  font "Times-Italic, 30"
set ylabel "P_2(P_{cm}) (GeV^{-3})" font "Times-Italic, 30"
set logscale y
set yr[1e1:3e4]
set xr[0: 1.2]
set ytics 10, 10, 1e4


set style line 1 lt 1 lw 5 lc rgb "red"
set style line 2 lt 1 lw 5 lc rgb "light-red"
set style line 3 lt 1 lw 5 lc rgb "navy"
set style line 4 lt 1 lw 5 lc rgb "green"
set style line 5 lt 1 lw 5 lc rgb "cyan"
set style line 6 lt 1 lw 5 lc rgb "web-green"


plot "hocommom.C" u 1:($2)/($1)/($1) ti "^{12}C all pairs" ls 3 w l, \
"datahofitall.C" u 1:($2)/($1)/($1)  w l ls 5  ti "{/Symbol s}= 140 MeV", \
"hocommom.C" u 1:($3)/($1)/($1) ti "^{12}C | ^1S_0" ls 1 w l, \
"datahofit.C" u 1:($2)/($1)/($1):($3)/($1)/($1) with filledcurves ls 4  ti "{/Symbol s}= 161 +/- 4 MeV "



