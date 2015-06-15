#!/usr/bin/gnuplot
set term postscript enhanced color "Helvetica, 29"
set output "eff_mass.eps"

set border lw 4
set tics scale 1,0.5


set key height 1 left
set nokey


#set style data line
set style line 2 lt 2 lw 5 lc rgb "#E41A1C" 
set style line 3 lt 3 lw 5 lc rgb "#377eb8"
set style line 4 lt 4 lw 5 lc rgb "#4daf4a"
set style line 5 lt 5 lw 5 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"



#set multiplot layout 2, 2 columnsfirst upwards

#set lmargin 8
#set rmargin 0.5
#set bmargin 2.8
#set tmargin 0.7

set ylabel "effective mass ratio {/Symbol h}" offset 2,0
set xlabel "(2/A) (N_{pn}+N_{pp}+N_{nn})"  offset 1.9,0
#set xlabel "per nucleon SRC probability"

set yr[1:1.04]
set xr[2:10]
#set xtics 0,1,6
set ytics 1, 0.01, 1.04

#set logscale y
#set yr [-0.05: 0.01]
#set format y "10^{%L}"
#set ytics -0.04, 0.02, 0.0 
#set mytics 2

set label '^{4}He' at 3-0.7,1.0115+0.003 font "Helvetica,25"
set label '^{9}Be' at 3.617-0.2,1.0125-0.003 font "Helvetica,25"
set label '^{12}C' at 4.277-0.5,1.0160+0.003 font "Helvetica,25"
set label '^{56}Fe' at 7.19-1.2,1.0255+0.002 font "Helvetica,25"
set label '^{197}Au' at 9.576-1.4,1.0330+0.002 font "Helvetica,25"

plot "eff_mass" u 4:2 ls 2 ps 2.5, 1.0012+ x*0.0033 ls 3


#0.07857588 + x* 0.08854221 ls 4

#unset multiplot


