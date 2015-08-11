#!/usr/bin/gnuplot
set term postscript eps enhanced color "Helvetica, 26"
set output "p2n1_p.eps"

set style data points

#set size 0.6,0.6
#set xtics 10
#set mxtics 10

set lmargin 9
set rmargin 2
set key left bottom

set format y "10^{%L}"
set logscale  y
set yr[1e-5:2e3]
set ytic 1e-5, 100, 1e3
set style data line


#set style line 3 lt 3 lw 5 lc rgb "#8c510a" 
#set style line 2 lt 2 lw 5 lc rgb "#d8b365" 
#set style line 1 lt 1 lw 5 lc rgb "#f6e8c3"
#set style line 4 lt 4 lw 5 lc rgb "#c7eae5"
#set style line 5 lt 5 lw 5 lc rgb "#5ab4ac"
#set style line 6 lt 6 lw 5 lc rgb "#01665e"

set xr[0:3]

set multiplot
#set key left top
set size 1,0.5
set origin 0.0,0.0

set bmargin 3
set tmargin 0

#set yr [0:1600]
#set ytic 300, 300, 1200
#set mytics 2

set xlabel "k (fm)" offset 0,0.5
set ylabel "HO k^2 n_1^{(0)}(k)"

plot "../result/ho_n1_p_0_3.C" u 1:($1)*($1)*($2) ti "C N" lw 5,\
"../result/ho_n1_p_0_3.Fe" u 1:($1)*($1)*($2) ti "Fe N" lw 5,\
"../result/ho_n1_p_0_3.Pb" u 1:($1)*($1)*($2) ti "Pb N" lw 5,\
"../result/ho_n1_p_0_3.C" u 1:($1)*($1)*($3) ti "C P" lw 5,\
"../result/ho_n1_p_0_3.Fe" u 1:($1)*($1)*($3) ti "Fe P" lw 5,\
"../result/ho_n1_p_0_3.Pb" u 1:($1)*($1)*($3) ti "Pb P" lw 5

set origin 0.0,0.5
set size 1,0.5

set bmargin 0
set tmargin 3

#set yr [0:1600]
set ytic 1e-3, 100, 1e3
#set mytics 10

set format x ""

set xlabel ""
set ylabel "WS k^2n_1^{(0)}(k)"

plot "../result/ws_n1_p_0_3.C" u 1:($1)*($1)*($2) ti "C N", "" u 1:($1)*($1)*($3) ti "C P", "../result/ws_n1_p_0_3.Fe" u 1:($1)*($1)*($2) ti "Fe N", "" u 1:($1)*($1)*($3) ti "Fe P", "../result/ws_n1_p_0_3.Pb" u 1:($1)*($1)*($2) ti "Pb N", "" u 1:($1)*($1)*($3) ti "Pb P"

set nomultiplot
