#!/usr/bin/bash

gnuplot << ??
set term epslatex standalone color  ",bx,14" dl 2.5 
set output "dens_rel_all.tex"

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

set bmargin 3.2
set lmargin 7
set tmargin 0.6
set rmargin 1.5
set style data points
set border lw 3
set tics scale 2,1

hbarc = 0.197327

set format y "\$\\\bm{10^{%L}}\$"
set xr[0:1]
set mytics 5
set yr[1e-1:1e7]
#set ytics 1e-3, 10, 5e3
set logscale y
set xlabel "\$\\\bm{k_{12}}\$ [GeV]"
set nokey


set ylabel "\$\\\bm{n_2(k_{12})}\$ [GeV\$^{\\\bm{-3}}\$]" offset 2.1, 0
set key samplen 1.5
set key top right height 1



plot "dens_rel.all.Fe" u (\$1)*hbarc:(\$2)/hbarc/hbarc/hbarc ls 2 w l ti 'mf' ,\
"" u (\$1)*hbarc:(\$3)/hbarc/hbarc/hbarc ls 3 w l ti 'cen' ,\
"" u (\$1)*hbarc:(\$4)/hbarc/hbarc/hbarc ls 4 w l ti 'ten' ,\
"" u (\$1)*hbarc:(\$5)/hbarc/hbarc/hbarc ls 5 w l ti 'iso' ,\
"" u (\$1)*hbarc:(\$6)/hbarc/hbarc/hbarc ls 6 w l ti 'all' 

??

pdflatex dens_rel_all.tex
rm dens_rel_all-inc.eps
rm dens_rel_all-inc-eps-converted-to.pdf
rm dens_rel_all.{aux,log}
