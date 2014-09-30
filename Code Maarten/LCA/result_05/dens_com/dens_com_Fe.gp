#!/usr/bin/bash
#
# ./dens_com.gp 
#

fig1="dens_com_$1_v1.eps"
fig2="dens_com_$1_v2.eps"

gnuplot << ??

set term postscript enhanced color "Computer-Modern-Roman, 30"

set style line 2 lt 1 lw 6 lc rgb "#E41A1C"
set style line 3 lt 1 lw 6 lc rgb "#377eb8"
set style line 4 lt 1 lw 6 lc rgb "#4daf4a"
set style line 5 lt 5 lw 6 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"

set bmargin 3.2
set lmargin 7 
set tmargin 0.35
set rmargin 1.2
set style data points
set border lw 2
set tics scale 1.5,0.75
#set key font ",25"

hbarc = 0.197327

set xr[0:0.4]

set xlabel "P_{12} [GeV]"
set ylabel "P_2(P_{12}|l=0) [GeV^{-3}]" offset 2.3

set grid front
unset grid



set output "$fig1"
plot "dens_com.$1.-10" u (\$1)*hbarc:-1:(\$2)/hbarc/hbarc/hbarc ls 2 w filledcu ti "0s-0s" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$2)/hbarc/hbarc/hbarc:(\$2+\$3)/hbarc/hbarc/hbarc ls 3 w filledcu ti "0s-0p" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$2+\$3)/hbarc/hbarc/hbarc:(\$2+\$3+\$4)/hbarc/hbarc/hbarc ls 4 w filledcu ti "0s-0d" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$5)/hbarc/hbarc/hbarc ls -1 lw 3  w l ti "total"

set output "$fig2"
plot "dens_com.$1.-10" u (\$1)*hbarc:(\$2)/hbarc/hbarc/hbarc ls 2 w l ti "0s-0s" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$3)/hbarc/hbarc/hbarc ls 3 w l ti "0s-0p" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$4)/hbarc/hbarc/hbarc ls 4 w l ti "0s-0d" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$5)/hbarc/hbarc/hbarc ls 5 w l ti "0s-0f" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$6)/hbarc/hbarc/hbarc ls 6 w l ti "0s-1s" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$7)/hbarc/hbarc/hbarc ls 7 w l ti "0p-0p" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$8)/hbarc/hbarc/hbarc ls 8 w l ti "0p-0d" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$9)/hbarc/hbarc/hbarc ls 9 w l ti "0p-0f" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$10)/hbarc/hbarc/hbarc ls 10 w l ti "0p-1s" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$11)/hbarc/hbarc/hbarc ls 11 w l ti "0d-0d" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$12)/hbarc/hbarc/hbarc ls 12 w l ti "0d-0f" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$13+\$15)/hbarc/hbarc/hbarc ls 13 w l ti "0d-1s" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$14)/hbarc/hbarc/hbarc ls 14 w l ti "0f-0f" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$16)/hbarc/hbarc/hbarc ls 15 w l ti "1s-0f" ,\
"dens_com.$1.-10" u (\$1)*hbarc:(\$17)/hbarc/hbarc/hbarc ls 16 w l ti "1s-1s" 

??
