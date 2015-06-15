#!/usr/bin/bash
# no ARG

file="dens_corr_ob_He"

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
set format y "\$\\\bm{10^{%L}}\$"
set mytics 5
set yr[1e-4:1e3]
set ytics 1e-4, 10, 1e3
set logscale y
set xlabel "\$\\\bm{k}\$ [fm\$\\\bm{^{-1}}\$]"
set nokey



set ylabel "\$\\\bm{n_1(k)}\$ " offset 2.1, 0
set key samplen 1.5
set key top right height 1 Right spacing 1.2



plot "dens_ob3.00.110.He.-1" u 1:(\$2+\$3) w l ls 2 ti '\$\bm{^{4}}\$He \$\bm{n^{[1]}(\vec{p})}\$' ,\
"dens_ob3.00.110.He.-1" u 1:(\$2) w l ls 3 ti 'IPM' ,\
"dens_ob3.00.100.He.-1" u 1:(\$3) w l  ls 4 ti 'central' ,\
"dens_ob3.00.010.He.-1" u 1:(\$3) w l  ls 5 ti 'tensor'


??

pdflatex $file
rm $file-inc-eps-converted-to.pdf
#rm $file-inc.eps
rm $file.{aux,log}
