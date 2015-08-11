#!/usr/bin/gnuplot

set term postscript enhanced color "Helvetica, 25"
set output "n2mf_D.eps"

set style data points
set border lw 2
set tics scale 1,0.5
set key font ",25"

set xr[0:6]
set yr[1e-4:1e2]
set xlabel "k [fm^{-1}]"
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

set ylabel "n_2^{mf}" offset 3, 0
plot "rel_dens_mf.D" ls 2 w l 

set output "n22b_D.eps"
set ylabel "n_2" offset 3, 0
set key
plot "rel_dens_mf.D" ls 2 w l ti "mf", "rel_dens_2b.D" ls 3 w l ti "2b" ,\
"rel_dens_2b_cen.D" ls 4 w l ti "2b central",\
"rel_dens_2b_ten.D" ls 5 w l ti "2b tensor",\
"rel_dens.D" ls 6 w l ti "total" ,\
"deuteron_p.dat"  u (($1)*sqrt(2)):(($4)/sqrt(8)) ls 7 ti "Real. D" w l
