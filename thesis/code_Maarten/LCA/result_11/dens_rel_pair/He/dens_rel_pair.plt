#!/usr/bin/bash

file="dens_rel_pair_$1"

echo $file

gnuplot << ??
set term epslatex standalone color  ",bx,17" dl 2.5 header '\usepackage{bm}'
set output "$file.tex"

set style line 1 lt 1 lw 4 lc rgb "#000000"
set style line 2 lt 2 lw 6 lc rgb "#E41A1C"
set style line 3 lt 3 lw 6 lc rgb "#377eb8"
set style line 4 lt 4 lw 6 lc rgb "#4daf4a"
set style line 5 lt 5 lw 6 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"
set style line 10 lt 2 lw 4 lc rgb "#000000"

set bmargin 2.6
set tmargin 0.1
set lmargin 6.9
set rmargin 0.6
#set style data points
set border lw 3
set tics scale 2,1

set format y "\$\\\bm{10^{%L}}\$"
set xr[0:5]
set yr[1e-6:5e1]
set mytics 5
set ytics 1e-6, 10, 5e1
set logscale y
set xlabel '\$\bm{k_{12}}\$ [fm\$\bm{^{-1}}\$]' offset 0,0.5


set ylabel "\$\\\bm{n_2(k_{12})}\$" offset 2.1, 0
set key samplen 1.5
set key top right height 1


norm=$2*($2-1)/2

plot "dens_rel.00.111.$1.-1" u (\$1):(\$4/norm) ls 2 w l ti '\$\bm{^{$2}}\$$1 all' ,\
"dens_rel.-11.111.$1.-1" u (\$1):(\$4/norm) ls 3 w l ti 'pn' ,\
"dens_rel.11.111.$1.-1" u (\$1):(\$4/norm) ls 4 w l ti 'pp' ,\
"dens_rel.-1-1.111.$1.-1" u (\$1):(\$4/norm) ls 5 w l ti 'nn' ,\
'$1$2_qNN.table.txt' u 1:((\$2)/norm/19.7392088) ls 7 ti 'Argonne pn',\
'$1$2_qNN.table.txt' u 1:((\$4)/norm/19.7392088) ls 9 ti 'Argonne pp/nn'

??

pdflatex $file
rm $file-inc.eps
rm $file-inc-eps-converted-to.pdf
rm $file.{aux,log}

