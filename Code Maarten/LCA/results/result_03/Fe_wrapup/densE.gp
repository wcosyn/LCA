#!/usr/bin/bash
#
# ./densE.gp X logymin logymax E 
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


set ylabel "n_2(k_{rel})" offset 3, 0
set key

#
# OUTDATED NEED TO BE UPDATED
#
#set output "momdens.all.$1.eps"
#
#plot "dens_mf.E-1.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 2 w l ti "mf",\
#"dens_2b.E-1.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 3 w l ti "2b" ,\
#"dens_2b_cen.E-1.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 4 w l ti "2b central",\
#"dens_2b_ten.E-1.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 5 w l ti "2b tensor",\
#"dens.E-1.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 6 w l ti "total" 
#
#set output "momdens.E$4.$1.eps"
#
#plot "dens_mf.E$4.$1" u ((\$1)/sqrt(2)):((\$4)*sqrt(8)) ls 2 w l ti "mf E$4",\
#"dens_2b.E$4.$1" u ((\$1)/sqrt(2)):((\$4)*sqrt(8)) ls 3 w l ti "2b E$4" ,\
#"dens_2b_cen.E$4.$1" u ((\$1)/sqrt(2)):((\$4)*sqrt(8)) ls 4 w l ti "2b central E$4",\
#"dens_2b_ten.E$4.$1" u ((\$1)/sqrt(2)):((\$4)*sqrt(8)) ls 5 w l ti "2b tensor E$4",\
#"dens.E$4.$1" u ((\$1)/sqrt(2)):((\$4)*sqrt(8)) ls 6 w l ti "total E$4" 
#

set output "dens_rel_all_E$4.$1.eps"

plot "dens_rel_all.E-1.$1" u 1:2 ls 2 w l ti "total" ,\
"dens_rel_all.E-1.$1" u 1:3 ls 3 w l ti "2b" ,\
"dens_rel_all.E$4.$1" u 1:3 ls 4 w l ti "2b E$4"

#"dens_2b_cen.E$4.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 4 w l ti "2b central E$4",\
#"dens_2b_ten.E$4.$1" u ((\$1)/sqrt(2)):((\$2)*sqrt(8)) ls 5 w l ti "2b tensor E$4"


??
