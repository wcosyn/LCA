#!/usr/bin/bash
# no ARG

file="corrf"

echo $file

gnuplot << ??

mpl_top    = 0.2 #inch  outer top margin, title goes here
mpl_bot    = 0.75 #inch  outer bottom margin, x label goes here
mpl_left   = 1.0 #inch  outer left margin, y label goes here
mpl_right  = 0.1 #inch  outer right margin, y2 label goes here
mpl_height = 4.0 #inch  height of individual plots
mpl_width  = 3.0 #inch  width of individual plots
mpl_dx     = 1.2 #inch  inter-plot horizontal spacing
mpl_dy     = 0.3 #inch  inter-plot vertical spacing
mpl_ny     = 1   #number of rows
mpl_nx     = 2   #number of columns

# calculate full dimensions
xsize = mpl_left+mpl_right+(mpl_width*mpl_nx)+(mpl_nx-1)*mpl_dx
ysize = mpl_top+mpl_bot+(mpl_ny*mpl_height)+(mpl_ny-1)*mpl_dy

# placement functions
#   rows are numbered from bottom to top
bot(n) = (mpl_bot+(n-1)*mpl_height+(n-1)*mpl_dy)/ysize
top(n)  = 1-((mpl_top+(mpl_ny-n)*(mpl_height+mpl_dy))/ysize)
#   columns are numbered from left to right
left(n) = (mpl_left+(n-1)*mpl_width+(n-1)*mpl_dx)/xsize
right(n)  = 1-((mpl_right+(mpl_nx-n)*(mpl_width+mpl_dx))/xsize)

set term epslatex standalone color  ",bx," 17 dl 5.1  size xsize,ysize header '\usepackage{bm}'
set output "$file".".tex"

hbarc = 0.197327

set style line 1 lt 1 lw 6 lc rgb "#000000"
set style line 2 lt 2 lw 6 lc rgb "#E41A1C"
set style line 3 lt 6 lw 6 lc rgb "#377eb8"
#set style line 4 lt 4 lw 6 lc rgb "#4daf4a"
set style line 4 lt 4 lw 6 lc rgb "#006600"
set style line 5 lt 1 lw 6 lc rgb "#984EA3"
set style line 6 lt 3 lw 6 lc rgb "#ff7f00"
set style line 7 lt 3 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 1 lw 6 lc rgb "#999999"
set style line 10 lt 2 lw 4 lc rgb "#000000"

set offsets
set autoscale fix
set size 1,1
set style data points
set border lw 4
set tics scale 2,1

# define axis setting for all subplots
set xr[0:3.5]
set xtics 0,1,5 offset -0.3
set mytics 5
set yr[-0.4:1.2]
set key samplen 2.2
set key top right height 1 Right spacing 1.2

kF=1.25


#start plotting
set multiplot

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set ylabel 'correlation function' offset 1, 0
#set format y "%1.1f"
set xlabel '\$\bm{r_{12} = \left|\vec{r}_1 - \vec{r}_2\right|}\$ [fm]'
#set format x ""



plot "./central_r.out" u (\$1):5 w l ls 2 ti 'central',\
  "./tensor_r.out" u (\$1):((\$3)*5) w l ls 3 ti 'tensor (x5)',\
  "./spinisospin_r.out" u (\$1):(\$2*5) w l ls 4 ti 'spin-isospin (x5)'

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set yr[1e-6:2e-2]
set xr[0:4]
set xtics 0,1,5 offset -0.3
set logscale y

set ylabel '\$\bm{ | f(p_{12}) |^2}\$ [fm\$\bm{^6}\$]' offset 1, 0
set format y "\$\\\bm{10^{%L}}\$"
set xlabel '\$\bm{p_{12}=\left|\frac{\vec{p}_1 - \vec{p}_2}{2}\right|}\$ [fm\$\bm{^{-1}}\$]'

set arrow from 1.25, 10**-6 to 1.25, 2e-2 nohead ls 9
plot "./central_p.out" u (\$1):((\$5)**2) w l ls 2 noti,\
  "./tensor_p.out" u (\$1):((\$3)**2) w l ls 3 noti,\
  "./spinisospin_p.out" u (\$1):((\$2)**2) w l ls 4 noti 


??

pdflatex $file
rm $file-inc-eps-converted-to.pdf
#rm $file-inc.eps
rm $file.{aux,log}
