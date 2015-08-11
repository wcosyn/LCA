#!/usr/bin/gnuplot

set term postscript enhanced color "Helvetica, 25"
set output "dens_mf.eps"

set style data points
set border lw 2
set tics scale 1,0.5
set key font ",25"

set xr[0:3]
set yr[1e-2:1e5]
set xlabel "k_{rel} [fm^{-1}]"
set logscale y

set style line 2 lt 2 lw 6 lc rgb "#E41A1C" 
set style line 3 lt 3 lw 6 lc rgb "#377eb8"
set style line 4 lt 4 lw 6 lc rgb "#4daf4a"
set style line 5 lt 5 lw 6 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"

set ylabel "n_2^{mf}(k_{rel})" offset 3, 0
plot "rel_dens_mf.D" u ($1)/sqrt(2):($2)*sqrt(8) ls 2 w l ti "D" ,\
"rel_dens_mf.He" u ($1)/sqrt(2):($2)*sqrt(8) ls 3 w l ti "He" ,\
"rel_dens_mf.C" u ($1)/sqrt(2):($2)*sqrt(8) ls 4 w l ti "C" ,\
"rel_dens_mf.Fe" u ($1)/sqrt(2):($2)*sqrt(8) ls 5 w l ti "Fe"

