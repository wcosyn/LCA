#!/usr/bin/bash
#
#set term postscript size 19,10.5 enhanced color "Helvetica, 25" 
gnuplot << ??
set output "moments_all.tex"
set term epslatex standalone color ",bx,11" header '\usepackage{bm}'

set style data points
set border lw 3
set tics scale 1,0.5
set key samplen 1
#set key font ",25"


set style line 2 lt 1 lw 4 ps 1.5 lc rgb "#E41A1C" 
set style line 3 lt 3 lw 4 ps 1.5 lc rgb "#377eb8"
set style line 4 lt 4 lw 4 ps 1.5 lc rgb "#4daf4a"
set style line 5 lt 5 lw 4 ps 1.5 lc rgb "#984EA3"
set style line 6 lt 6 lw 4 ps 1.5 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"



#set offsets graph 0, 0, 0.005, 0
set multiplot layout 2, 2 columnsfirst downwards 


#set key top right
set key top right
set label 4 "\$\\\bm{(a)\\\; ^{12}}\$C" at -0.75, 0.195 left

set ylabel "width [GeV]" offset 2, -5

set style data histogram
set style fill solid noborder
set style histogram clustered

set yr[0.09:0.21]
set ytics 0.1, 0.02, 0.2
set xr[-1:6]
unset xtics
set ytics nooffset
set lmargin 7
set rmargin 0
set bmargin 0
set tmargin 3
plot for [COL=2:5] 'momentsC.x.f2' u COL:xticlabels(1) title columnheader ls COL,\
"results_tang.dat" i 0 u (0.6):(\$2)/1000:(\$3)/1000 with yerrorbars noti lt -1 pt 2 ps 1.7 ,\
"wsmomentsC.x" i 0 u (\$1)+1.35:2 with points ti "WS pp" ls 6,\
"" i 1 u ((\$1)+1.35):2 with points noti ls 6

#"WSmoments.dat" i 0 u (1.35):(\$2)/1000 with points ti "WS pp" ls 6
#"results_Gfits.dat" i 0 u (1.35):(\$4)/1000:(\$5)/1000 with yerrorbars lt -1 pt 3 ps 1.7 ti "WS pp (G)"
unset label 4

set label 1 "\$\\\bm{(b)\\\; ^{27}}\$Al" at -0.75, 0.195 left
set nokey #ti  "^{27}Al"
unset ylabel
set xtics nomirror
set lmargin 7
set rmargin 0
set tmargin 0
set bmargin 3
set xtics rotate 
plot for [COL=2:5] 'momentsAl.x.f2' u COL:xticlabels(1) ls COL title columnheader,\
"wsmomentsAl.x" i 0 u (\$1)+1.35:2 with points ti "WS pp" ls 6,\
"" i 1 u ((\$1)+1.35):2 with points noti ls 6

#"WSmoments.dat" i 1 u (1.35):(\$2)/1000 with points  ls 6 noti,\
#"results_Gfits.dat" i 1 u (1.35):(\$4)/1000:(\$5)/1000 with yerrorbars lt -1 pt 3 ps 1.7 ti "G pp"
unset label 1

set label 2 "\$\\\bm{(c)\\\; ^{56}}\$Fe" at -0.75, 0.195 left
set xr[-1:8]
#set key ti  "^{56}Fe"
set ytics ("" 0.1, "" 0.12, "" 0.14, "" 0.16, "" 0.18, "" 0.2)
unset xtics
set rmargin 7
set lmargin 0
set bmargin 0
set tmargin 3
plot for [COL=2:5] 'momentsFe.x.f2' u COL:xticlabels(1) ls COL title columnheader,\
"wsmomentsFe.x" i 0 u (\$1)+1.35:2 with points ti "WS pp" ls 6,\
"" i 1 u ((\$1)+1.35):2 with points noti ls 6

#"WSmoments.dat" i 2 u (1.35):(\$2)/1000 with points ls 6 noti,\
#"results_Gfits.dat" i 2 u (1.35):(\$4)/1000:(\$5)/1000 with yerrorbars lt -1 pt 3 ps 1.7 ti "G pp"
unset label 2

set label 3 "\$\\\bm{(d)\\\; ^{208}}\$Pb" at -0.75, 0.195 left
set nokey
#set key ti  "^{208}Pb"
set xtics nomirror rotate
set rmargin 7
set lmargin 0
set tmargin 0
set bmargin 3
plot for [COL=2:5] 'momentsPb.x.f2' u COL:xticlabels(1) ls COL title columnheader,\
"wsmomentsPb.x" i 0 u (\$1)+1.35:2 with points ti "WS pp" ls 6,\
"" i 1 u ((\$1)+1.35):2 with points noti ls 6

#"WSmoments.dat" i 3 u (1.35):(\$2)/1000 with points  ls 6 noti,\
#"results_Gfits.dat" i 3 u (1.35):(\$4)/1000:(\$5)/1000 with yerrorbars lt -1 pt 3 ps 1.7 ti "G pp"
unset multiplot
??
pdflatex moments_all.tex


