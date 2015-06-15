#!/usr/bin/bash
#
# ./densKin.gp X
#

gnuplot << ??
set term postscript enhanced color "Helvetica, 25"

set style line 2 lt 2 lw 6 lc rgb "#E41A1C" 
set style line 3 lt 3 lw 6 lc rgb "#377eb8"
set style line 4 lt 4 lw 6 lc rgb "#4daf4a"
set style line 5 lt 5 lw 6 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"

set style data points
set border lw 2
set tics scale 1,0.5
set key font ",25"

set xr[0:3.5]
set xlabel "k_{rel} [fm^{-1}]"


set ylabel "n_2(k_{rel})" offset 3, 0
set key


set output "dens_rel_kin_all.$1.eps"

plot "dens_rel_kin_all.E-1.$1" u 1:(\$3)/(\$2) ls 2 w l ti "mf" ,\
"dens_rel_kin_all.E-1.$1" u 1:(\$4)/(\$2) ls 3 w l ti "2b" ,\
"dens_rel_kin_all.E0.$1" u 1:5 ls 4 w l ti "2b E0",\
"dens_rel_kin_all.E1.$1" u 1:5 ls 5 w l ti "2b E1",\
"dens_rel_kin_all.E2.$1" u 1:5 ls 6 w l ti "2b E2",\
"dens_rel_kin_all.E3.$1" u 1:5 ls 7 w l ti "2b E3"

??
