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

#set style data points
set border lw 2
set tics scale 1,0.5
set key font ",24"

#set yr[0:4]
set xr[0:4]
#set logscale y
#set key bottom right


set xlabel "k_{rel} [fm^{-1}]"
#set ylabel "{/Symbol f}_{tensor}(k_{rel})"
set ylabel "R_t^{(nl)}"

set output "overlap_tensor.all.$1.eps"
plot "overlap.00.$1" u 1:3 ls 2 w l ti "n=0,l=0" ,\
 "overlap.all.$1" u 1:3 ls 3 w l ti "all" 

set output "overlap_tensor.l0.$1.eps"
plot "overlap.00.$1" u 1:3 ls 2 w l ti "n=0,l=0" ,\
"overlap.10.$1" u 1:3 ls 3 w l ti "n=1,l=0" ,\
"overlap.20.$1" u 1:3 ls 4 w l ti "n=2,l=0" ,\
"overlap.30.$1" u 1:3 ls 5 w l ti "n=3,l=0" ,\
"overlap.40.$1" u 1:3 ls 6 w l ti "n=4,l=0"

set output "overlap_tensor.l1.$1.eps"
plot "overlap.01.$1" u 1:3 ls 2 w l ti "n=0,l=1" ,\
"overlap.11.$1" u 1:3 ls 3 w l ti "n=1,l=1" ,\
"overlap.21.$1" u 1:3 ls 4 w l ti "n=2,l=1" ,\
"overlap.31.$1" u 1:3 ls 5 w l ti "n=3,l=1" ,\
"overlap.41.$1" u 1:3 ls 6 w l ti "n=4,l=1"

set output "overlap_tensor.l2.$1.eps"
plot "overlap.02.$1" u 1:3 ls 2 w l ti "n=0,l=2" ,\
"overlap.12.$1" u 1:3 ls 3 w l ti "n=1,l=2" ,\
"overlap.22.$1" u 1:3 ls 4 w l ti "n=2,l=2" ,\
"overlap.32.$1" u 1:3 ls 5 w l ti "n=3,l=2" ,\
"overlap.42.$1" u 1:3 ls 6 w l ti "n=4,l=2"

set output "overlap_tensor.l3.$1.eps"
plot "overlap.03.$1" u 1:3 ls 2 w l ti "n=0,l=3" ,\
"overlap.13.$1" u 1:3 ls 3 w l ti "n=1,l=3" ,\
"overlap.23.$1" u 1:3 ls 4 w l ti "n=2,l=3" ,\
"overlap.33.$1" u 1:3 ls 5 w l ti "n=3,l=3" ,\
"overlap.43.$1" u 1:3 ls 6 w l ti "n=4,l=3"

#set ylabel "{/Symbol f}_{central}(k_{rel})"
set ylabel "R_t^{(nl)}"

set output "overlap_central.all.$1.eps"
plot "overlap.00.$1" u 1:2 ls 2 w l ti "n=0,l=0" ,\
 "overlap.all.$1" u 1:2 ls 3 w l ti "all" 

set output "overlap_central.l0.$1.eps"
plot "overlap.00.$1" u 1:2 ls 2 w l ti "n=0,l=0" ,\
"overlap.10.$1" u 1:2 ls 3 w l ti "n=1,l=0" ,\
"overlap.20.$1" u 1:2 ls 4 w l ti "n=2,l=0" ,\
"overlap.30.$1" u 1:2 ls 5 w l ti "n=3,l=0" ,\
"overlap.40.$1" u 1:2 ls 6 w l ti "n=4,l=0"

set output "overlap_central.l1.$1.eps"
plot "overlap.01.$1" u 1:2 ls 2 w l ti "n=0,l=1" ,\
"overlap.11.$1" u 1:2 ls 3 w l ti "n=1,l=1" ,\
"overlap.21.$1" u 1:2 ls 4 w l ti "n=2,l=1" ,\
"overlap.31.$1" u 1:2 ls 5 w l ti "n=3,l=1" ,\
"overlap.41.$1" u 1:2 ls 6 w l ti "n=4,l=1"

set output "overlap_central.l2.$1.eps"
plot "overlap.02.$1" u 1:2 ls 2 w l ti "n=0,l=2" ,\
"overlap.12.$1" u 1:2 ls 3 w l ti "n=1,l=2" ,\
"overlap.22.$1" u 1:2 ls 4 w l ti "n=2,l=2" ,\
"overlap.32.$1" u 1:2 ls 5 w l ti "n=3,l=2" ,\
"overlap.42.$1" u 1:2 ls 6 w l ti "n=4,l=2"

set output "overlap_central.l3.$1.eps"
plot "overlap.03.$1" u 1:2 ls 2 w l ti "n=0,l=3" ,\
"overlap.13.$1" u 1:2 ls 3 w l ti "n=1,l=3" ,\
"overlap.23.$1" u 1:2 ls 4 w l ti "n=2,l=3" ,\
"overlap.33.$1" u 1:2 ls 5 w l ti "n=3,l=3" ,\
"overlap.43.$1" u 1:2 ls 6 w l ti "n=4,l=3"

??
