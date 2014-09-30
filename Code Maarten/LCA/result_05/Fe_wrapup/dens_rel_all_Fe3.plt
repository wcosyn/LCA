#!/usr/bin/bash
gnuplot << ??
set term epslatex standalone color  ",bx,14" dl 2.5 
set output "dens_rel_all_Fe3.tex"
#set term postscript enhanced color "Helvetica, 25"
#set term postscript enhanced color "Computer-Modern-Roman, 30"

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
#set key font ",25"

hbarc = 0.197327

#set xr[0:3.5*hbarc]
set format y "\$\\\bm{10^{%L}}\$"
set xr[0:0.7]
set mytics 5
set yr[1e-3:2e3]
set ytics 1e-3, 10, 5e3
set logscale y
set xlabel "\$\\\bm{k_{12}}\$ [GeV]"
set nokey


#set ylabel "n_2^{2n+l}(k) / n_2^{all}(k)" offset 1.5, 0
set ylabel "\$\\\bm{n_2(k_{12})}\$ [GeV\$^{\\\bm{-3}}\$]" offset 2.1, 0
set key samplen 1.5
set key top right height 1



plot "dens_rel_all.E0.Fe3" u (\$1)*hbarc:(\$2)/hbarc/hbarc/hbarc/1558 ls 2 w l ti "\$\\\bm{2n+l = 0}\$" ,\
"dens_rel_all.E1.Fe3" u (\$1)*hbarc:(\$2)/hbarc/hbarc/hbarc/1558 ls 3 w l ti "\$\\\bm{2n+l = 1}\$" ,\
"dens_rel_all.E2+.Fe3" u (\$1)*hbarc:(\$2)/hbarc/hbarc/hbarc/1558 ls 4 w l ti "\$\\\bm{2n+l \\\geq 2}\$" ,\
"dens_rel_all.E-1.Fe3" u (\$1)*hbarc:(\$2)/hbarc/hbarc/hbarc/1558 ls 1 w l ti "TBC",\
"dens_rel_all.E-1.Fe3" u (\$1)*hbarc:(\$2-\$3)/hbarc/hbarc/hbarc/1558 ls 10 w l ti "IPM"

#"dens_rel_all.E3.Fe3" u (\$1)*hbarc:2 ls 5 w l ti "2n+l = 3" ,\
#"dens_rel_all.E4.Fe3" u (\$1)*hbarc:2 ls 6 w l ti "2n+l = 4" ,\
#"dens_rel_all.E5.Fe3" u (\$1)*hbarc:2 ls 7 w l ti "2n+l = 5" ,\
#"dens_rel_all.E6.Fe3" u (\$1)*hbarc:2 ls 8 w l ti "2n+l = 6" ,\

??

pdflatex dens_rel_all_Fe3.tex
