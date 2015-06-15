#!/usr/bin/gnuplot
set term postscript size 15,10.5 enhanced color "Helvetica, 30" 
set output "moments_pp_vs_A.eps"

set style data points
set border lw 2
set tics scale 1,0.5
set key samplen 2
set key top center


set style line 2 lt 2 lw 6 ps 1.5 lc rgb "#E41A1C" 
set style line 3 lt 3 lw 6 ps 1.5 lc rgb "#377eb8"
set style line 4 lt 4 lw 6 ps 1.5 lc rgb "#4daf4a"
set style line 5 lt 5 lw 6 ps 1.5 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 ps 1.5 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"

set logscale x
set xr [10:250]
set xlabel "A"
set ylabel "{/Symbol s}_{cm} [GeV]"




plot "moments_pp_vs_A.dat" u 1:2 w l ls 2 ti "l=0",\
"" u 1:3 w l ls 3 ti "l=1",\
"" u 1:4 w l ls 4 ti  "l=2",\
"" u 1:($5)/1000:($6)/1000 with errorbars pt 2 lt -1 ps 3 ti "experiment"
