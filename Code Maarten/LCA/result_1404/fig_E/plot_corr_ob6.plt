#!/usr/bin/bash
# no ARG

file="dens_ob_E"

echo $file

gnuplot << ??

mpl_top    = 0.1 #inch  outer top margin, title goes here
mpl_bot    = 0.75 #inch  outer bottom margin, x label goes here
mpl_left   = 1.0 #inch  outer left margin, y label goes here
mpl_right  = 0.1 #inch  outer right margin, y2 label goes here
mpl_height = 2.4 #inch  height of individual plots
mpl_width  = 3.2 #inch  width of individual plots
mpl_dx     = 0.3 #inch  inter-plot horizontal spacing
mpl_dy     = 0.3 #inch  inter-plot vertical spacing
mpl_ny     = 2   #number of rows
mpl_nx     = 3   #number of columns

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

set term epslatex standalone color  ",bx," 17 dl 4  size xsize,ysize header '\usepackage{bm}'
set output "$file".".tex"

hbarc = 0.197327

#set style line 1 lt 1 lw 6 lc rgb "#000000"
set style line 2 lt 2 lw 6 lc rgb "#E41A1C"
set style line 3 lt 6 lw 6 lc rgb "#377eb8"
#set style line 4 lt 1 lw 6 lc rgb "#4daf4a"
set style line 4 lt 1 lw 6 lc rgb "#006600"
set style line 5 lt 2 lw 6 lc rgb "#984EA3"
set style line 6 lt 4 lw 6 lc rgb "#ff7f00"
set style line 7 lt 3 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#999999"
set style line 9 lt 9 lw 6 lc rgb "#999999"
set style line 10 lt 2 lw 4 lc rgb "#000000"

set offsets
set autoscale fix
set size 1,1
set style data points
set border lw 4
set tics scale 2,1

# define axis setting for all subplots
set xr[0:5]
set xtics 0,1,5 offset -0.3
set format y "\$\\\bm{10^{%L}}\$"
set mytics 5
set yr[1e-4:2e2]
set ytics 1e-4, 10, 1e3 offset 0, -0.25
set logscale y
set key samplen 1.9
set key top right height 1 Right spacing 1.05

set pointsize 1.2

f(x) = (x<0) ? 1e-6:x

#start plotting
set multiplot

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(2)
set bmargin at screen bot(2)

set ylabel '\$\bm{n^{[1]}(p)}\$ [fm\$\bm{^3}\$]' offset 2.1, 0
set xlabel ""
set format x ""

A=4
X="He"

#plot "he4.momentum.txt" i 0 u 1:((\$2)*2/19.7392088/A) ls 5 ti '\$\bm{^{'.A.'}}\$'.X.' Argonne' ,\

#plot "he4.momentum.txt" i 0 u 1:((\$2)*2)/19.7392088/A  ls 5 noti,\


plot "./dens_ob4.00.110.".X.".-1-1-1-1" u 1:(f(\$4)/A) w l ls 2 ti '\$\bm{^{'.A.'}}\$'.X.' LCA' ,\
"" u 1:(f(\$2)/A) w l  ls 3 ti 'IPM',\
"./dens_ob4.00.110.".X.".0000" u 1:(f(\$3)/A) w l  ls 4 ti '\$\bm{00,00}\$'


#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(2)
set bmargin at screen bot(2)

set ylabel ""
set format y ""

A=9
X="Be"

#plot "c12.momentum.txt" i 0 u 1:((\$2)*2/19.7392088/A/1.10)  ls 5 noti ,\

plot "./dens_ob4.00.110.".X.".-1-1-1-1" u 1:(f(\$4)/A) w l ls 2 ti '\$\bm{^{'.A.'}}\$'.X.' LCA' ,\
"" u 1:(f(\$2)/A) w l  ls 3 ti 'IPM',\
"./dens_ob4.00.110.".X.".0000" u 1:(f(\$3)/A) w l  ls 4 ti '\$\bm{00,00}\$',\
"./dens_ob4.00.110.".X.".diff" u 1:(f(\$3)/A) w l  ls 5 ti 'Other'

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(3)
set rmargin at screen right(3)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(2)
set bmargin at screen bot(2)

set ylabel ""
set format y ""

A=16
X="O"


plot "./dens_ob4.00.110.".X.".-1-1-1-1" u 1:(f(\$4)/A) w l ls 2 ti '\$\bm{^{'.A.'}}\$'.X.' LCA' ,\
"" u 1:(f(\$2)/A) w l  ls 3 ti 'IPM',\
"./dens_ob4.00.110.".X.".0000" u 1:(f(\$3)/A) w l  ls 4 ti '\$\bm{00,00}\$',\
"./dens_ob4.00.110.".X.".diff" u 1:(f(\$3)/A) w l  ls 5 ti 'Other'

#plot "../".X."/dens_ob3.00.110.".X.".-1" u 1:(f(\$4)/A) w l ls 2 ti '\$\bm{^{'.A.'}}\$'.X.' LCA' ,\
#  "../".X."/dens_ob3.00.110.".X.".-1" u 1:(f(\$2)/A) w l  ls 3 ti 'IPM',\
#"../".X."/dens_ob3.00.110.".X.".0" u 1:(f(\$3)/A) w l  ls 4 ti '\$\bm{2n+l=0}\$',\
#"../".X."/dens_ob3.00.110.".X.".1" u 1:(f(\$3)/A) w l  ls 6 ti '\$\bm{2n+l=1}\$',\
#"../".X."/dens_ob3.00.110.".X.".2" u 1:(f(\$3)/A) w l  ls 7 ti '\$\bm{2n+l=2}\$',\
#"../".X."/dens_ob3.00.110.".X.".3+" u 1:(f(\$3)/A) w l  ls 8 ti '\$\bm{2n+l\geq3}\$'

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set ylabel '\$\bm{n^{[1]}(p)}\$ [fm\$\bm{^3}\$]' offset 2.1, 0
set format y "\$\\\bm{10^{%L}}\$"
set xlabel "\$\\\bm{p}\$ [fm\$\\\bm{^{-1}}\$]"
set format x "%1.0f"

A=27
X="Al"
Xkey="Al"


plot "./dens_ob4.00.110.".X.".-1-1-1-1" u 1:(f(\$4)/A) w l ls 2 ti '\$\bm{^{'.A.'}}\$'.X.' LCA' ,\
"" u 1:(f(\$2)/A) w l  ls 3 ti 'IPM',\
"./dens_ob4.00.110.".X.".0000" u 1:(f(\$3)/A) w l  ls 4 ti '\$\bm{00,00}\$',\
"./dens_ob4.00.110.".X.".1010" u 1:(f(\$3)/A) w l  ls 6 ti '\$\bm{10,10}\$',\
"./dens_ob4.00.110.".X.".diff" u 1:(f(\$3)/A) w l  ls 5 ti 'Other'

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set ylabel ""
set format y ""
set xlabel "\$\\\bm{p}\$ [fm\$\\\bm{^{-1}}\$]"
set format x "%1.0f"

A=48
X="Ca48"
Xkey="Ca"

plot "./dens_ob4.00.110.".X.".-1-1-1-1" u 1:(f(\$4)/A) w l ls 2 ti '\$\bm{^{'.A.'}}\$'.Xkey.' LCA' ,\
"" u 1:(f(\$2)/A) w l  ls 3 ti 'IPM',\
"./dens_ob4.00.110.".X.".0000" u 1:(f(\$3)/A) w l  ls 4 ti '\$\bm{00,00}\$',\
"./dens_ob4.00.110.".X.".1010" u 1:(f(\$3)/A) w l  ls 6 ti '\$\bm{10,10}\$',\
"./dens_ob4.00.110.".X.".diff" u 1:(f(\$3)/A) w l  ls 5 ti 'Other'

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(3)
set rmargin at screen right(3)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set ylabel ""
set format y ""

A=108
X="Ag"


plot "../".X."/dens_ob4.00.110.".X.".-1" u 1:(f(\$4)/A) w l ls 2 ti '\$\bm{^{'.A.'}}\$'.X.' LCA' ,\
  "../".X."/dens_ob4.00.110.".X.".-1" u 1:(f(\$2)/A) w l  ls 3 ti 'IPM',\
"./dens_ob4.00.110.".X.".0000" u 1:(f(\$3)/A) w l  ls 4 ti '\$\bm{00,00}\$',\
"./dens_ob4.00.110.".X.".1010" u 1:(f(\$3)/A) w l  ls 6 ti '\$\bm{10,10}\$',\
"./dens_ob4.00.110.".X.".2020" u 1:(f(\$3)/A) w l  ls 7 ti '\$\bm{20,20}\$',\
"./dens_ob4.00.110.".X.".diff" u 1:(f(\$3)/A) w l  ls 5 ti 'Other'



??

pdflatex $file
rm $file-inc-eps-converted-to.pdf
#rm $file-inc.eps
rm $file.{aux,log}
