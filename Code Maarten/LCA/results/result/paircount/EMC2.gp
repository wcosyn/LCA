#!/usr/bin/gnuplot
set term postscript enhanced color "Helvetica, 29"
set output "EMC2.eps"

set border lw 4
set tics scale 1,0.5


set key height 1 left
set nokey


#set style data line
set style line 2 lt 2 lw 5 lc rgb "#E41A1C" 
set style line 3 lt 3 lw 5 lc rgb "#377eb8"
set style line 4 lt 4 lw 5 lc rgb "#4daf4a"
set style line 5 lt 5 lw 5 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"



#set multiplot layout 2, 2 columnsfirst upwards

#set lmargin 8
#set rmargin 0.5
#set bmargin 2.8
#set tmargin 0.7

set ylabel "-dR/dx_B" offset 0,0.4
set xlabel "(2/A) (N_{pn} + N_{pp} + N_{nn})"  offset 1.9,0
#set xlabel "per nucleon SRC probability"
load "EMClabels2.gp"

set yr[0.12:0.6]
set xr[2:10]
#set xtics 0,1,6

#set logscale y
#set yr [-0.05: 0.01]
#set format y "10^{%L}"
#set ytics -0.04, 0.02, 0.0 
#set mytics 2

plot "EMC" u ($6+$7+$8+$9)*2/($1):3:4 with errorbars pt 12 ps 2 lw 4,\
 0.11396 + x* 0.036446 ls 3

#plot "EMC" u ($6+$7+$8+$9)*4/($1)/($1-1):3:4 with errorbars pt 12 ps 2 lw 4



#0.07857588 + x* 0.08854221 ls 4

#unset multiplot


