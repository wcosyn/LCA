#!/usr/bin/bash

file="dens_ang_Cpp"

echo $file

gnuplot << ??
set term epslatex standalone color  ",bx,14" dl 2.5 
set output "$file".".tex"

#set term postscript enhanced color "Helvetica, 25"
#set term postscript enhanced color "Computer-Modern-Roman, 30"

set style line 2 lt 2 lw 6 lc rgb "#E41A1C"
set style line 3 lt 3 lw 6 lc rgb "#377eb8"
set style line 4 lt 4 lw 6 lc rgb "#4daf4a"

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
set yr[1e-3:1e3]
#set ytics 1e-3, 10, 1e6
set logscale y
set xlabel "\$\\\bm{k_{12}}\$ [fm\$\\\bm{^{-1}}\$]"
set nokey



set ylabel '\$\bm{n_2(k_{12},Q, \cos \theta )}\$' offset 2.1, 0
set key samplen 1.5
set key top right height 1 spacing 2
set label '\$\bm{Q=0}\$ fm\$\bm{^{-1}}\$' at 0.25, 3e2 left
set label '\$\bm{Q=1}\$ fm\$\bm{^{-1}}\$' at 2, 1e-2 left
set label '\$\bm{Q=2}\$ fm\$\bm{^{-1}}\$' at 0.25, 1e-2 left




plot "dens_ang.Cpp.0" i 0 u 2:5 w l ls 2 ti '\$\bm{\cos \theta =1}\$' ,\
"" i 20 u 2:5 w l ls 2  noti,\
"" i 40 u 2:5 w l ls 2  noti,\
"dens_ang.Cpp.45" i 0 u 2:5 w l ls 3 ti '\$\bm{\cos\theta=\frac{\sqrt{2}}{2}}\$',\
"" i 20 u 2:5 w l ls 3  noti,\
"" i 40 u 2:5 w l ls 3  noti,\
"dens_ang.Cpp.90" i 0 u 2:5 w l ls 4 ti '\$\bm{\cos\theta=0}\$' ,\
"" i 20 u 2:5 w l ls 4  noti,\
"" i 40 u 2:5 w l ls 4  noti 

??

pdflatex $file
rm $file-inc-eps-converted-to.pdf
rm $file-inc.eps
rm $file.{aux,log}
