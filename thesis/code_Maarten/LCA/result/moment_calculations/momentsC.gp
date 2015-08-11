#!/usr/bin/gnuplot

set term postscript enhanced color "Helvetica, 20"
set output "momentsC.eps"

set style data points
set border lw 2
set tics scale 1,0.5
set key font ",18"

set xtics 0,1,2
set xr[-1.5:2.5]
set xtics ("All" -1, "l = 0" 0, "1" 1, "2" 2)
set ylabel "width (GeV)"


set style line 2 lt 2 lw 5 lc rgb "#E41A1C" 
set style line 3 lt 3 lw 5 lc rgb "#377eb8"
set style line 4 lt 4 lw 5 lc rgb "#4daf4a"
set style line 5 lt 5 lw 5 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"

set multiplot layout 2, 2

plot "momentsC.x" i 0 u ($1)-0.05:2 ti "nn" ls 2 ,\
"" i 1 u ($1)-0.05:2 noti ls 2 ,\
"" i 2 u ($1)-0.025:2 ti "np S=0" ls 3 ,\
"" i 3 u ($1)-0.025:2 noti ls 3 ,\
"" i 4 u ($1)+0.025:2 ti "np S=1" ls 4 ,\
"" i 5 u ($1)+0.025:2 noti ls 4 ,\
"" i 6 u ($1)+0.05:2 ti "pp" ls 5,\
"" i 7 u ($1)+0.05:2 noti ls 5,\
"experiment.dat" i 0 u (-0.5):($2)/1000:($3)/1000 with yerrorbars lt -1  noti

plot "momentsC.y" i 0 u ($1)-0.05:2 ti "nn" ls 2 ,\
"" i 1 u ($1)-0.05:2 noti ls 2 ,\
"" i 2 u ($1)-0.025:2 ti "np S=0" ls 3 ,\
"" i 3 u ($1)-0.025:2 noti ls 3 ,\
"" i 4 u ($1)+0.025:2 ti "np S=1" ls 4 ,\
"" i 5 u ($1)+0.025:2 noti ls 4 ,\
"" i 6 u ($1)+0.05:2 ti "pp" ls 5,\
"" i 7 u ($1)+0.05:2 noti ls 5,\
"experiment.dat" i 0 u (-0.5):($2)/1000:($3)/1000 with yerrorbars lt -1  noti

plot "momentsC.y" i 0 u ($1)-0.05:2 ti "nn" ls 2 ,\
"" i 1 u ($1)-0.05:2 noti ls 2 ,\
"" i 2 u ($1)-0.025:2 ti "np S=0" ls 3 ,\
"" i 3 u ($1)-0.025:2 noti ls 3 ,\
"" i 4 u ($1)+0.025:2 ti "np S=1" ls 4 ,\
"" i 5 u ($1)+0.025:2 noti ls 4 ,\
"" i 6 u ($1)+0.05:2 ti "pp" ls 5,\
"" i 7 u ($1)+0.05:2 noti ls 5,\
"experiment.dat" i 0 u (-0.5):($2)/1000:($3)/1000 with yerrorbars lt -1  noti

unset multiplot
