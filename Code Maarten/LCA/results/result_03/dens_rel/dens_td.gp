#!/usr/bin/bash
#
# ./dens_td.gp X logymin logymax E 
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
set logscale y
set yr[1e$2:1e$3]
set nokey


set ylabel "n_2(k_{rel}, P_{cm}=0)" offset 3, 0
set key

set output "dens_td_all_E$4.$1.eps"

plot "dens_td_all.E-1.$1" i 0 u 2:3 ls 2 w l ti "total" ,\
"dens_td_all.E-1.$1" i 0 u 2:4 ls 3 w l ti "2b" ,\
"dens_td_all.E$4.$1" i 0 u 2:4 ls 4 w l ti "2b E$4"

#"dens_2b_cen.E$4.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 4 w l ti "2b central E$4",\
#"dens_2b_ten.E$4.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 5 w l ti "2b tensor E$4"


??
