#!/bin/bash

gnuplot << ??
set term postscript enhanced color "Helvetica, 25"

set style data points
set border lw 2
set tics scale 1,0.5
set key font ",25"

set xr[0:4]
#set yr[1e-4:1e2]
set yr[1e$2:1e$3]
set xlabel "k_{rel} [fm^{-1}]"
set logscale y
set nokey

set style line 2 lt 2 lw 6 lc rgb "#E41A1C" 
set style line 3 lt 3 lw 6 lc rgb "#377eb8"
set style line 4 lt 4 lw 6 lc rgb "#4daf4a"
set style line 5 lt 5 lw 6 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"

set output "n22b_$1.eps"
set ylabel "n_2(k_{rel})" offset 3, 0
set key
plot "rel_dens_mf.E-1.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 2 w l ti "mf",\
"rel_dens_2b.E-1.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 3 w l ti "2b" ,\
"rel_dens_2b_cen.E-1.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 4 w l ti "2b central",\
"rel_dens_2b_ten.E-1.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 5 w l ti "2b tensor",\
"rel_dens.E-1.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 6 w l ti "total" 
#"deuteron_p.dat"  u 1:4 ls 7 ti "Real. D" w l
#"deuteron_p.dat"  u ((\$1)*sqrt(2)):((\$4)/sqrt(8)) ls 7 ti "Real. D" w l

??
