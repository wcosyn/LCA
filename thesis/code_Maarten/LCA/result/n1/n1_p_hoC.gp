#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color "Helvetica, 29"
set output "n1_p_hoC.eps"

set style data line
set key left top height 1
set border lw 2
set tics scale 1,1
set pointsize 2.5
set key spacing 2

set xr[0:320]
set xlabel "k [fm]" offset 0, 0.2
set ylabel "k^2 n^{(0)}_1(k)" offset 2, 0
set xtics 0, 50, 300

set style line 1 lt 1 lw 5 lc rgb "#8c510a"
set style line 2 lt 1 lw 5 lc rgb "#d8b365"
set style line 3 lt 1 lw 5 lc rgb "#b2182b"
set style line 4 lt 1 lw 5 lc rgb "#ef8a62"

plot \
"hofits.C" u ($1)*197.327:2:3 with filledcurve ls 2 noti, \
"../result/ho_n1_p_0_3.C" u ($1)*197.327:($1)*($1)*($4) ls 1 ti "^{12}C HO"

