#!/usr/bin/gnuplot
#set term postscript size 15,15.75 enhanced color "Helvetica, 25" 

c=5 
t=6

plot "overlap.27" u 2:( ($1)==0 ? (column(c)) : 1/0 ) ti "n=0" ls 2,\
"" u 2:( ($1)==1 ? column(c) : 1/0 ) ti "n=1" ls 3,\
"" u 2:( ($1)==2 ? column(c) : 1/0 ) ti "n=2" ls 4,\
"" u 2:( ($1)==3 ? column(c) : 1/0 ) ti "n=3" ls 5,\
"" u 2:( ($1)==4 ? column(c) : 1/0 ) ti "n=4" ls 6,\
"" u 2:( ($1)==5 ? column(c) : 1/0 ) ti "n=5" ls 7

MAX1=GPVAL_Y_MAX
MIN1=GPVAL_Y_MIN

plot "overlap.27" u 2:( ($1)==0 ? column(t) : 1/0 ) ti "n=0" ls 2,\
"" u 2:( ($1)==1 ? column(t) : 1/0 ) ti "n=1" ls 3,\
"" u 2:( ($1)==2 ? column(t) : 1/0 ) ti "n=2" ls 4,\
"" u 2:( ($1)==3 ? column(t) : 1/0 ) ti "n=3" ls 5,\
"" u 2:( ($1)==4 ? column(t) : 1/0 ) ti "n=4" ls 6,\
"" u 2:( ($1)==5 ? column(t) : 1/0 ) ti "n=5" ls 7

MAX2=GPVAL_Y_MAX
MIN2=GPVAL_Y_MIN

set term postscript size 15,5.25 enhanced color "Helvetica, 25" 
set output "overlaps.eps"

set style data points
set border lw 2
set tics scale 1,0.5
set key samplen 1
set key center right
set key font ",25"


set style line 2 lt 2 lw 4 ps 1.5 lc rgb "#E41A1C" 
set style line 3 lt 3 lw 4 ps 1.5 lc rgb "#377eb8"
set style line 4 lt 4 lw 4 ps 1.5 lc rgb "#4daf4a"
set style line 5 lt 5 lw 4 ps 1.5 lc rgb "#984EA3"
set style line 6 lt 6 lw 4 ps 1.5 lc rgb "#ff7f00"
set style line 7 lt 7 lw 4 ps 1.5 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"


#set multiplot layout 3, 2 rowsfirst downwards
set multiplot layout 1, 2 rowsfirst downwards

#set lmargin 8
#set rmargin 1
#set bmargin 4
#set tmargin 1


set xr[-0.5:5.5]
set xtics ( "l = 0" 0, "1" 1, "2" 2, "3" 3, "4" 4, "5" 5)

set ylabel "Al Central Corr Strength" offset 2, 0
#set yr[-0.025:0]
#set ytics -0.02, 0.01, 0

set yrange [MIN1-(MAX1-MIN1)*0.10:MAX1+(MAX1-MIN1)*0.10]

plot "overlap.27" u 2:( ($1)==0 ? (column(c)) : 1/0 ) ti "n=0" ls 2,\
"" u 2:( ($1)==1 ? column(c) : 1/0 ) ti "n=1" ls 3,\
"" u 2:( ($1)==2 ? column(c) : 1/0 ) ti "n=2" ls 4,\
"" u 2:( ($1)==3 ? column(c) : 1/0 ) ti "n=3" ls 5,\
"" u 2:( ($1)==4 ? column(c) : 1/0 ) ti "n=4" ls 6,\
"" u 2:( ($1)==5 ? column(c) : 1/0 ) ti "n=5" ls 7


set ylabel "Al Tensor Corr Strength" offset 2, 0
#set yr[-0.045:0]
#set ytics -0.04, 0.01, 0

set key top right
set yrange [MIN2-(MAX2-MIN2)*0.10:MAX2+(MAX2-MIN2)*0.10]

plot "overlap.27" u 2:( ($1)==0 ? column(t) : 1/0 ) ti "n=0" ls 2,\
"" u 2:( ($1)==1 ? column(t) : 1/0 ) ti "n=1" ls 3,\
"" u 2:( ($1)==2 ? column(t) : 1/0 ) ti "n=2" ls 4,\
"" u 2:( ($1)==3 ? column(t) : 1/0 ) ti "n=3" ls 5,\
"" u 2:( ($1)==4 ? column(t) : 1/0 ) ti "n=4" ls 6,\
"" u 2:( ($1)==5 ? column(t) : 1/0 ) ti "n=5" ls 7

exit

set ylabel "Fe Central Corr Strength" offset 2, 0
#set yr[-0.025:0]
#set ytics -0.02, 0.01, 0

plot "overlap.56" u 2:( ($1)==0 ? column(c) : 1/0 ) ti "n=0" ls 2,\
"" u 2:( ($1)==1 ? column(c) : 1/0 ) ti "n=1" ls 3,\
"" u 2:( ($1)==2 ? column(c) : 1/0 ) ti "n=2" ls 4,\
"" u 2:( ($1)==3 ? column(c) : 1/0 ) ti "n=3" ls 5,\
"" u 2:( ($1)==4 ? column(c) : 1/0 ) ti "n=4" ls 6,\
"" u 2:( ($1)==5 ? column(c) : 1/0 ) ti "n=5" ls 7

set ylabel "Fe Tensor Corr Strength" offset 2, 0
#set yr[-0.04:0]
#set ytics -0.04, 0.01, 0

plot "overlap.56" u 2:( ($1)==0 ? column(t) : 1/0 ) ti "n=0" ls 2,\
"" u 2:( ($1)==1 ? column(t) : 1/0 ) ti "n=1" ls 3,\
"" u 2:( ($1)==2 ? column(t) : 1/0 ) ti "n=2" ls 4,\
"" u 2:( ($1)==3 ? column(t) : 1/0 ) ti "n=3" ls 5,\
"" u 2:( ($1)==4 ? column(t) : 1/0 ) ti "n=4" ls 6,\
"" u 2:( ($1)==5 ? column(t) : 1/0 ) ti "n=5" ls 7

set ylabel "Pb Central Corr Strength" offset 2, 0
#set yr[-0.02:0]
#set ytics -0.02, 0.01, 0

plot "overlap.208" u 2:( ($1)==0 ? column(c) : 1/0 ) ti "n=0" ls 2,\
"" u 2:( ($1)==1 ? column(c) : 1/0 ) ti "n=1" ls 3,\
"" u 2:( ($1)==2 ? column(c) : 1/0 ) ti "n=2" ls 4,\
"" u 2:( ($1)==3 ? column(c) : 1/0 ) ti "n=3" ls 5,\
"" u 2:( ($1)==4 ? column(c) : 1/0 ) ti "n=4" ls 6,\
"" u 2:( ($1)==5 ? column(c) : 1/0 ) ti "n=5" ls 7

set ylabel "Pb Tensor Corr Strength" offset 2, 0
#set yr[-0.035:0]
#set ytics -0.04, 0.01, 0

plot "overlap.208" u 2:( ($1)==0 ? column(t) : 1/0 ) ti "n=0" ls 2,\
"" u 2:( ($1)==1 ? column(t) : 1/0 ) ti "n=1" ls 3,\
"" u 2:( ($1)==2 ? column(t) : 1/0 ) ti "n=2" ls 4,\
"" u 2:( ($1)==3 ? column(t) : 1/0 ) ti "n=3" ls 5,\
"" u 2:( ($1)==4 ? column(t) : 1/0 ) ti "n=4" ls 6,\
"" u 2:( ($1)==5 ? column(t) : 1/0 ) ti "n=5" ls 7
