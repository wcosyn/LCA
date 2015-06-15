#!/usr/bin/gnuplot
set term postscript size 15,10.5 enhanced color "Helvetica, 25" 
set output "moments_all_l0.eps"

set style data points
set border lw 2
set tics scale 1,0.5
set key samplen 1
#set key font ",25"


set style line 2 lt 2 lw 4 ps 1.5 lc rgb "#E41A1C" 
set style line 3 lt 3 lw 4 ps 1.5 lc rgb "#377eb8"
set style line 4 lt 4 lw 4 ps 1.5 lc rgb "#4daf4a"
set style line 5 lt 5 lw 4 ps 1.5 lc rgb "#984EA3"
set style line 6 lt 6 lw 4 ps 1.5 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"


set multiplot layout 2, 2 rowsfirst 

set lmargin 8
set rmargin 1
set bmargin 2
set tmargin 1

set xtics 0,1,2
set xr[-1.5:1.5]
set xtics ("l=0" -1, "n = 0" 0, "1" 1)
#set xlabel "C"
set ylabel "width (GeV)" offset 2, 0
set label  "^{12}C" at graph 0.5,0.2 
#set title "C"

plot "momentsC.x" i 0 u ($1)-0.1:2 ti "nn" ls 2 ,\
"" i 1 u ($1)-0.1:2 noti ls 2 ,\
"" i 2 u ($1)-0.05:2 ti "np S=0" ls 3 ,\
"" i 3 u ($1)-0.05:2 noti ls 3 ,\
"" i 4 u ($1)+0.05:2 ti "np S=1" ls 4 ,\
"" i 5 u ($1)+0.05:2 noti ls 4 ,\
"" i 6 u ($1)+0.1:2 ti "pp" ls 5,\
"" i 7 u ($1)+0.1:2 noti ls 5,\
"WSmoments.dat" i 0 u (0):($2)/1000 ti "WS pp" ls 6,\
"experiment.dat" i 0 u (-0.5):($2)/1000:($3)/1000 with yerrorbars noti lt -1 pt 10 ps 1.7 ,\
"results_tang.dat" i 0 u (-0.5):($2)/1000:($3)/1000 with yerrorbars noti lt -1 pt 2 ps 1.7 
# "results_Gfits.dat" i 0 u 0:($2)/1000:($3)/1000 with yerrorbars lt -1 pt 3 ps 1.7 ti "G pp"

set nokey

set xtics 0,1,2
set xr[-1.5:2.5]
set xtics ("l=0" -1, "n = 0" 0, "1" 1, "2" 2)
#set xlabel "Al"
set ylabel "width (GeV)"
unset label
set label  "^{27}Al" at graph 0.3,0.2 
#set title "Al"

plot "momentsAl.x" i 0 u ($1)-0.1:2 ti "nn" ls 2 ,\
"" i 1 u ($1)-0.1:2 noti ls 2 ,\
"" i 2 u ($1)-0.05:2 ti "np S=0" ls 3 ,\
"" i 3 u ($1)-0.05:2 noti ls 3 ,\
"" i 4 u ($1)+0.05:2 ti "np S=1" ls 4 ,\
"" i 5 u ($1)+0.05:2 noti ls 4 ,\
"" i 6 u ($1)+0.1:2 ti "pp" ls 5,\
"" i 7 u ($1)+0.1:2 noti ls 5,\
"WSmoments.dat" i 1 u 0.0:($2)/1000 ti "WS pp" ls 6,\
"experiment.dat" i 1 u (-0.5):($2)/1000:($3)/1000 with yerrorbars noti lt -1 pt 10 ps 1.7
#"results_Gfits.dat" i 1 u 0:($2)/1000:($3)/1000 with yerrorbars lt -1 pt 3 ps 1.7 ti "G pp"


set xtics 0,1,2
set xr[-1.5:3.5]
set xtics ("l=0" -1, "n = 0" 0, "1" 1, "2" 2, "3" 3 )
#set xlabel "Fe"
set ylabel "width (GeV)"
unset label
set label  "^{56}Fe" at graph 0.3,0.2 
#set title "Fe"
set nokey

plot "momentsFe.x" i 0 u ($1)-0.1:2 ti "nn" ls 2 ,\
"" i 1 u ($1)-0.1:2 noti ls 2 ,\
"" i 2 u ($1)-0.05:2 ti "np S=0" ls 3 ,\
"" i 3 u ($1)-0.05:2 noti ls 3 ,\
"" i 4 u ($1)+0.05:2 ti "np S=1" ls 4 ,\
"" i 5 u ($1)+0.05:2 noti ls 4 ,\
"" i 6 u ($1)+0.1:2 ti "pp" ls 5,\
"" i 7 u ($1)+0.1:2 noti ls 5,\
"WSmoments.dat" i 2 u 0.0:($2)/1000 ti "WS pp" ls 6,\
"experiment.dat" i 2 u (-0.5):($2)/1000:($3)/1000 with yerrorbars noti lt -1 pt 10 ps 1.7
# "results_Gfits.dat" i 2 u 0:($2)/1000:($3)/1000 with yerrorbars lt -1 pt 3 ps 1.7 ti "G pp"


set xtics 0,1,2
set xr[-1.5:4.5]
set xtics ("l=0" -1, "n = 0" 0, "1" 1, "2" 2, "3" 3, "4" 4 )
#set xlabel "Pb"
set ylabel "width (GeV)"
unset label
set label  "^{208}Pb" at graph 0.3,0.2 
#set title "Pb"
set nokey

plot "momentsPb.x" i 0 u ($1)-0.1:2 ti "nn" ls 2 ,\
"" i 1 u ($1)-0.1:2 noti ls 2 ,\
"" i 2 u ($1)-0.05:2 ti "np S=0" ls 3 ,\
"" i 3 u ($1)-0.05:2 noti ls 3 ,\
"" i 4 u ($1)+0.05:2 ti "np S=1" ls 4 ,\
"" i 5 u ($1)+0.05:2 noti ls 4 ,\
"" i 6 u ($1)+0.1:2 ti "pp" ls 5,\
"" i 7 u ($1)+0.1:2 noti ls 5,\
"WSmoments.dat" i 3 u 0.0:($2)/1000 ti "WS pp" ls 6,\
"experiment.dat" i 3 u (-0.5):($2)/1000:($3)/1000 with yerrorbars noti lt -1 pt 10 ps 1.7
# "results_Gfits.dat" i 3 u 0:($2)/1000:($3)/1000 with yerrorbars lt -1 pt 3 ps 1.7 ti "G pp"


unset multiplot


