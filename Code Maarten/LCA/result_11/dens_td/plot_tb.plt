#!/usr/bin/bash
# no ARG

file="dens_td_$1"
filemf="dens_td_mf_$1"



echo $file
echo $filemf

gnuplot << ??
set term epslatex standalone color ",bx," 17 dl 2.5 header '\usepackage{bm}'
set output "$file".".tex"

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

set border lw 3
set tics scale 2,1
set size ratio 1

set xr[0:5]
set yr[0:5]
set format cb "\$\\\bm{10^{%L}}\$"
#set mytics 5
#set yr[1e-4:1e3]
#set ytics 1e-4, 10, 1e3
#set logscale y
set xlabel "\$\\\bm{P}\$ [fm\$\\\bm{^{-1}}\$]"
set ylabel "\$\\\bm{k}\$ [fm\$\\\bm{^{-1}}\$]"
set cblabel "\$\\\bm{^{$2}}\$$1 ".'$\bm{n_2(k,P)}$'
set logscale cb
set cbr[1e-6:1e2]
set cbtics 1e-6,100,1e2
unset mcbtics

unset key
#set nokey
#set ylabel "\$\\\bm{n_1(k)}\$ " offset 2.1, 0
#set key samplen 1.5
#set key top right height 1 Right

set pm3d map
splot "dens_td.00.111.$1" u 1:2:(\$5/($2*($2-1))*2) 

set output "$filemf".".tex"

splot "dens_td.00.111.$1" u 1:2:(\$3/($2*($2-1))*2) 

??

pdflatex $file
pdflatex $filemf
rm $file-inc-eps-converted-to.pdf
rm $file-inc.eps
rm $file.{aux,log}
rm $filemf-inc-eps-converted-to.pdf
rm $filemf-inc.eps
rm $filemf.{aux,log}
