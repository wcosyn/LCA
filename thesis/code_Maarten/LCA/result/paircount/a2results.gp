#!/usr/bin/gnuplot
set term postscript eps color enhanced "Helvetica, 29"
set output "a2results.eps"

set border lw 4

set tics scale 2,1

set style line 2 lt 1 lw 5 lc rgb "#E41A1C" 
set style line 3 lt 1 lw 5 lc rgb "#377eb8"
set style line 4 lt 1 lw 4 lc rgb "#4daf4a"
set style line 5 lt 1 lw 4 lc rgb "#984EA3"
set style line 6 lt 1 lw 4 lc rgb "#ff7f00"
set style line 7 lt 1 lw 4 lc rgb "#a65628"
set style line 8 lt 1 lw 4 lc rgb "#ffff33"
set style line 9 lt 1 lw 4 lc rgb "#999999"

set style fill transparent solid 0.4 noborder


#set xrange[-20:8]
set log x
#set log y
set xr[1:250]
set yrange[0.9:8]
set xtics 1, 10, 300
set xlabel "mass number A" font "Times-Italic, 30"
set ylabel "a_{2} (A/d) " font "Times-Italic, 30"


set key right bottom height 1 
#set nokey

plot  "PairsNP.a2" u 1:4:5 w filledcurves ls 9 ti "with c.o.m",\
 "" u 1:3 ls 9 pt 13 ps 2 ti "without c.o.m." ,\
 "fomin" u 1:2:3 ti "Fomin et. al." with yerrorbars ls 2 ps 2 pt 2,\
 "egiyan" u 1:2:3 ti "Egiyan et. al." with yerrorbars ls 3 ps 2 pt 3,\
 "daya2.txt" u ($1):3:4 ti "SLAC" with yerrorbars ls 4 ps 2 pt 4,\
# "larrya2.txt" u 1:3:4 ti "JLAB (Hall B)" with yerrorbars ls 5 ps 2 pt 5,\
#  "newa2.txt" u ($1):3:4 ti "JLAB (Hall C)" with yerrorbars ls 6 ps 2 pt 6


