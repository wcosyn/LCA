#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color "Helvetica, 29"
set output "n1_p_wsFe.eps"

set style data line
set key left top height 1
set border lw 2
set tics scale 1,1
set pointsize 2.5
set key spacing 2

set xr[0:320]
set xlabel "k [MeV]" offset 0, 0.2
set ylabel "k^2 n^{(0)}_1(k)" offset 2, 0
set xtics 0, 50, 300

set style line 1 lt 1 lw 5 lc rgb "#8c510a"
set style line 2 lt 1 lw 5 lc rgb "#d8b365"
set style line 3 lt 1 lw 5 lc rgb "#b2182b"
set style line 4 lt 1 lw 5 lc rgb "#ef8a62"
set style line 5 lt 1 lw 5 lc rgb "#7fbf7b"
set style line 6 lt 1 lw 5 lc rgb "#1b7837"
hbarc= 197.327

plot \
"wsfits_kT.Fe" u ($1):($2):($3) with filledcurve ls 4 noti, \
"" u ($1):($2):($4) with filledcurve ls 4 noti, \
"" u ($1):($3):($4) with filledcurve ls 4 noti, \
"ws_n1_p_0_3.Fe" u ($1)*hbarc:($1)*($1)*($4)/hbarc ls 3 ti "^{56}Fe WS"

#plot \
#"wsfits.Fe" u ($1)*hbarc:($2)/hbarc:($3)/hbarc with filledcurve ls 4 noti, \
#"" u ($1)*hbarc:($2)/hbarc:($4)/hbarc with filledcurve ls 4 noti, \
#"" u ($1)*hbarc:($3)/hbarc:($4)/hbarc with filledcurve ls 4 noti, \
#"ws_n1_p_0_3.Fe" u ($1)*hbarc:($1)*($1)*($4)/hbarc ls 3 ti "^{56}Fe WS"
