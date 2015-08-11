#!/usr/bin/bash
# no ARG

file="dens_xx_ratio_$1"

echo $file

gnuplot << ??
set term epslatex standalone color  ",bx," 17 dl 2.5 header '\usepackage{bm}'
set output "$file".".tex"

#set term postscript enhanced color "Helvetica, 25"
#set term postscript enhanced color "Computer-Modern-Roman, 30"

set style line 1 lt 1 lw 4 lc rgb "#000000"
set style line 2 lt 2 lw 6 lc rgb "#E41A1C"
set style line 3 lt 4 lw 6 lc rgb "#377eb8"
set style line 4 lt 3 lw 6 lc rgb "#4daf4a"
set style line 5 lt 3 lw 6 lc rgb "#984EA3"
set style line 6 lt 3 lw 6 lc rgb "#ff7f00"
set style line 7 lt 3 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"
set style line 10 lt 2 lw 4 lc rgb "#000000"

set bmargin 3.2
set lmargin 7
set tmargin 0.6
set rmargin 1.5
set style data points
set border lw 3
set tics scale 2,1

hbarc = 0.197327

set xr[0:5]
#set format y "\$\\\bm{10^{%L}}\$"
set mytics 2
set yr[0:1]
set ytics 0, 0.2, 1
#set logscale y
set xlabel "\$\\\bm{p}\$ [fm\$\\\bm{^{-1}}\$]"



set ylabel '\$\bm{n^{[1]}_{NN}(p) \; / \; n^{[1]}(p)}\$' offset 2.1, 0
set key samplen 1.5
set key center right height 1 Right spacing 1.2



plot "dens_ratio.110.$1.-1" u 1:2 w l ls 2 ti 'pp',\
"dens_ratio.110.$1.-1" u 1:3 w l ls 3 ti 'nn',\
"dens_ratio.110.$1.-1" u 1:4 w l ls 4 ti 'pn' 


??

pdflatex $file
rm $file-inc-eps-converted-to.pdf
rm $file-inc.eps
rm $file.{aux,log}
