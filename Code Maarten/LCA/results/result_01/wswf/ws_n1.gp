#!/usr/bin/gnuplot
set term postscript size 15,5.25 enhanced color "Helvetica, 25" 
set output "ws_n1_C.eps"

set style data lines
set border lw 3
set tics scale 1,0.5
set key samplen 1
#set key font ",25"


set style line 2 lt 2 lw 6 ps 1.5 lc rgb "#E41A1C" 
set style line 3 lt 3 lw 6 ps 1.5 lc rgb "#377eb8"
set style line 4 lt 4 lw 6 ps 1.5 lc rgb "#4daf4a"
set style line 5 lt 5 lw 6 ps 1.5 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 ps 1.5 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"

set lmargin 8
set rmargin 1
set bmargin 4
set tmargin 1

set xr[0:2]
set xlabel "k [fm^{-1}]"

set multiplot layout 1, 2 rowsfirst 

set ylabel "n_2(k)" offset 2, 0
plot "ws_n1_p_0_5.C" i 0 u 1:($4)/12 ti "^{12}C" ls 2 

set nokey
set ylabel "k^2 n_2(k)" offset 2, 0
plot "ws_n1_p_0_5.C" i 0 u 1:(($1)*($1)*($4)/12) ls 2 

unset multiplot

set output "ws_n1_Al.eps"
set multiplot layout 1, 2 rowsfirst 

set key
set ylabel "n_2(k)" offset 2, 0
plot "ws_n1_p_0_5.Al" i 0 u 1:($4)/27 ti "^{27}Al" ls 5

set nokey
set ylabel "k^2 n_2(k)" offset 2, 0
plot "ws_n1_p_0_5.Al" i 0 u 1:(($1)*($1)*($4)/27) ls 5

unset multiplot

set output "ws_n1_Fe.eps"
set multiplot layout 1, 2 rowsfirst 

set key
set ylabel "n_2(k)" offset 2, 0
plot "ws_n1_p_0_5.Fe" i 0 u 1:($4)/56 ti "^{56}Fe" ls 3

set nokey
set ylabel "k^2 n_2(k)" offset 2, 0
plot "ws_n1_p_0_5.Fe" i 0 u 1:(($1)*($1)*($4)/56) ls 3

unset multiplot

set output "ws_n1_Pb.eps"
set multiplot layout 1, 2 rowsfirst 

set key
set ylabel "n_2(k)" offset 2, 0
plot "ws_n1_p_0_5.Pb" i 0 u 1:($4)/208 ti "^{208}Pb" ls 4

set nokey
set ylabel "k^2 n_2(k)" offset 2, 0
plot "ws_n1_p_0_5.Pb" i 0 u 1:(($1)*($1)*($4)/208) ls 4

unset multiplot
