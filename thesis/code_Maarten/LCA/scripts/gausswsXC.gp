#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color "Helvetica, 29"
set output "gaussXwsC.eps"

set style data line
set key right top height 1
set border lw 2
set tics scale 1,1
set pointsize 2.5
set key spacing 2

set style line 1 lt 1 lw 5 lc rgb "red"
set style line 2 lt 1 lw 5 lc rgb "light-red"
set style line 3 lt 1 lw 5 lc rgb "navy"
set style line 4 lt 1 lw 5 lc rgb "green"
set style line 5 lt 1 lw 5 lc rgb "cyan"
set style line 6 lt 1 lw 5 lc rgb "web-green"

set xlabel "P_{cm} (GeV)"  font "Times-Italic, 30"
set ylabel "P_2(P_{cm}) (GeV^{-2})" font "Times-Italic, 30"
set logscale y
set yr[1e1:2e5]
set xr[0:1]

plot "wscommom.C" u 1:($2)/($1)/($1) ti "^{12}C all pairs WS" ls 3 w l, \
"datahofitall.C" u 1:($2)/($1)/($1) ls 5  ti "{/Symbol s}= 140 MeV", \
"wscommom.C" u 1:($3)/($1)/($1) ti "^{12}C | ^1S_0 WS" ls 1 w l, \
"datawsfit.C" u 1:($2)/($1)/($1):($3)/($1)/($1) with filledcurves ls 4  ti "{/Symbol s}= 164 +/- 4 MeV"
