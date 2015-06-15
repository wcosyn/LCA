#!/usr/bin/bash

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
set yr[1e$2:1e$3]
set xlabel "k_{rel} [fm^{-1}]"
set logscale y
set nokey


set output "momdens.all.$1.eps"
set ylabel "n_2(k_{rel})" offset 3, 0
set key

plot "rel_dens_mf.n-1l-1.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 2 w l ti "mf",\
"rel_dens_2b.n-1l-1.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 3 w l ti "2b" ,\
"rel_dens_2b_cen.n-1l-1.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 4 w l ti "2b central",\
"rel_dens_2b_ten.n-1l-1.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 5 w l ti "2b tensor",\
"rel_dens.n-1l-1.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 6 w l ti "total" 

set output "momdens.n$4l$5.$1.eps"

plot "rel_dens_mf.n$4l$5.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 2 w l ti "mf n$4l$5",\
"rel_dens_2b.n$4l$5.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 3 w l ti "2b n$4l$5" ,\
"rel_dens_2b_cen.n$4l$5.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 4 w l ti "2b central n$4l$5",\
"rel_dens_2b_ten.n$4l$5.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 5 w l ti "2b tensor n$4l$5",\
"rel_dens.n$4l$5.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 6 w l ti "total n$4l$5" 

#unset logscale y
#set yr[0:2]
set output "momdens.n$4l$5.all.$1.eps"

plot "rel_dens.n-1l-1.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 2 w l ti "total" ,\
"rel_dens_2b.n-1l-1.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 3 w l ti "2b" ,\
"rel_dens_2b.n$4l$5.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 4 w l ti "2b n$4l$5",\
"rel_dens_2b_cen.n-1l-1.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 5 w l ti "2b central",\
"rel_dens_2b_cen.n$4l$5.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 6 w l ti "2b central n$4l$5",\
"rel_dens_2b_ten.n-1l-1.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 7 w l ti "2b tensor",\
"rel_dens_2b_ten.n$4l$5.$1" u ((\$1)/sqrt(2)):((\$3)*sqrt(8)) ls 9 w l ti "2b tensor n$4l$5"

??
