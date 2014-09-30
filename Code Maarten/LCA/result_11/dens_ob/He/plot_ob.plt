#!/usr/bin/bash

file="dens_ob_He"

echo $file

gnuplot << ??
set term epslatex standalone color  ",bx,17" dl 2.5 header '\usepackage{bm}'
set output "$file".".tex"

#set term postscript enhanced color "Helvetica, 25"
#set term postscript enhanced color "Computer-Modern-Roman, 30"

set style line 1 lt 1 lw 4 lc rgb "#000000"
set style line 2 lt 2 lw 6 lc rgb "#E41A1C"
set style line 3 lt 3 lw 6 lc rgb "#377eb8"
set style line 5 lt 2 lw 6 lc rgb "#4daf4a"
set style line 4 lt 3 lw 6 lc rgb "#984EA3"
set style line 6 lt 2 lw 6 lc rgb "#ff7f00"
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

set xr[0:5]
set format y "\$\\\bm{10^{%L}}\$"
set mytics 5
set yr[1e-4:1e3]
set ytics 1e-4, 10, 1e3
set logscale y
set xlabel "\$\\\bm{k}\$ [fm\$\\\bm{^{-1}}\$]"
set nokey



set ylabel "\$\\\bm{n_1(k)}\$ " offset 2.1, 0
set key samplen 1.5
set key top right height 1 Right



plot "he4.momentum.txt" u 1:(\$2/19.7392088*2) ls 2 ti "\$^4\$He Argonne - p+n" ,\
"dens_ob.He" u 1:(\$2+\$3) w l ls 3 ti "TBC - p+n" 

#"" u 1:4 ls 3 ti "pp" ,\
#"dens_td.-11.$3" i $5 u 2:(\$5*389.636) w l ls 5 ti "pn - TBC"
#"dens_td.pp.$3" i index u 2:(\$3*700) w l ls 4 ti "pp - MF" ,\
#"dens_td.pn.$3" i index u 2:(\$3*700) w l ls 6 ti "pn - MF" ,\

??

pdflatex $file
rm $file-inc-eps-converted-to.pdf
rm $file-inc.eps
rm $file.{aux,log}
