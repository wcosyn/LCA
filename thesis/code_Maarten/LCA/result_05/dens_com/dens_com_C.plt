#!/usr/bin/bash
#
# ./dens_com.gp 
#

#fig1="dens_com_$1.tex"
fig2="dens_com_C.tex"

gnuplot << ??

set term epslatex standalone color header '\usepackage{bm}' ',bx,14'  dl 3.5

#set term postscript enhanced color "Computer-Modern-Roman, 30"

set style line 2 lt 2 lw 6 lc rgb "#E41A1C"
set style line 3 lt 4 lw 6 lc rgb "#377eb8"
set style line 4 lt 5 lw 6 lc rgb "#4daf4a"
set style line 5 lt 5 lw 6 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"

set bmargin 3.2
set lmargin 7 
set tmargin 0.5
set rmargin 1.5
set style data points
set border lw 3
set tics scale 1.5,0.75
#set key font ",25"
set key spacing 1.5 
set key height 2 

hbarc = 0.197327

set xr[0:0.3]

set xlabel "\$\\\bm{P_{12}}\$ [GeV]"
set ylabel "\$\\\bm{P_2(P_{12}|nl=00)}\$ [GeV\$\\\bm{^{-3}}\$]" offset 1

set grid front
unset grid



#set output "$fig1"
#plot "dens_com.$1.-10" u (\$1)*hbarc:-1:(\$2)/hbarc/hbarc/hbarc ls 2 w filledcu ti "0s-0s" ,\
#"dens_com.$1.-10" u (\$1)*hbarc:(\$2)/hbarc/hbarc/hbarc:(\$2+\$3)/hbarc/hbarc/hbarc ls 3 w filledcu ti "0s-0p" ,\
#"dens_com.$1.-10" u (\$1)*hbarc:(\$2+\$3)/hbarc/hbarc/hbarc:(\$2+\$3+\$4)/hbarc/hbarc/hbarc ls 4 w filledcu ti "0p-0p" ,\
#"dens_com.$1.-10" u (\$1)*hbarc:(\$5)/hbarc/hbarc/hbarc ls -1 lw 3  w l ti "total"

set output "$fig2"
plot "dens_com.C.00" u (\$1)*hbarc:(\$2)/hbarc/hbarc/hbarc/66 ls 2 w l ti '$\bm{(0s0s)}$' ,\
  "dens_com.C.00" u (\$1)*hbarc:(\$3)/hbarc/hbarc/hbarc/66 ls 3 w l ti '$\bm{(0s0p)}$' ,\
  "dens_com.C.00" u (\$1)*hbarc:(\$4)/hbarc/hbarc/hbarc/66 ls 4 w l ti '$\bm{(0p0p)}$' ,\
"dens_com.C.00" u (\$1)*hbarc:(\$5)/hbarc/hbarc/hbarc/66 ls -1 lw 3  w l ti "total"
??

pdflatex dens_com_C.tex 
