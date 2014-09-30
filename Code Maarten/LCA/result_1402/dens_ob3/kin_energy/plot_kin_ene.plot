#!/usr/bin/bash
# no ARG

file="kin_cum"

echo $file

gnuplot << ??

mpl_top    = 0.1 #inch  outer top margin, title goes here
mpl_bot    = 0.75 #inch  outer bottom margin, x label goes here
mpl_left   = 1.0 #inch  outer left margin, y label goes here
mpl_right  = 0.1 #inch  outer right margin, y2 label goes here
mpl_height = 3.3 #inch  height of individual plots
mpl_width  = 4.4 #inch  width of individual plots
mpl_dx     = 0.3 #inch  inter-plot horizontal spacing
mpl_dy     = 0.3 #inch  inter-plot vertical spacing
mpl_ny     = 1   #number of rows
mpl_nx     = 1   #number of columns

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

set term epslatex standalone color  ",bx," 17 dl 3.5  size xsize,ysize header '\usepackage{bm}'
set output "$file".".tex"

hbarc = 0.197327

set style line 1 lt 1 lw 6 lc rgb "#000000"
set style line 2 lt 2 lw 6 lc rgb "#E41A1C"
set style line 3 lt 6 lw 6 lc rgb "#377eb8"
set style line 4 lt 4 lw 6 lc rgb "#4daf4a"
set style line 5 lt 1 lw 6 lc rgb "#984EA3"
set style line 6 lt 3 lw 6 lc rgb "#ff7f00"
set style line 7 lt 3 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"
set style line 10 lt 2 lw 4 lc rgb "#000000"

set offsets
set autoscale fix
set size 1,1
set style data points
set border lw 4
set tics scale 2,1

# define axis setting for all subplots
set xr[0:4]
set xtics 0,1,5 offset -0.3
set mytics 5
set yr[0:35]
set key samplen 0.5
set key top left height 1 Right spacing 1.2
set xlabel "\$\\\bm{p}\$ [fm\$\\\bm{^{-1}}\$]"



#start plotting
#set multiplot

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set ylabel 'cum. \$\bm{<T>}\$ [MeV]' offset 2.1, 0

A=12
X="C"

set label 1 at 2.5, 27 "TOTAL"
set label 3 at 2, 13.5 "IPM"
set label 2 at 2, 5 "CORR"

plot "kin_cum.12" u 1:(\$2) w l lt 2 lw 6   lc rgb "#E41A1C" ti '\$\bm{^{12}}\$C p/n',\
"" u 1:(\$3) w l lt 3 lw 6  lc rgb "#E41A1C" noti,\
"" u 1:(\$4) w l lt 4 lw 6  lc rgb "#E41A1C" noti,\
"kin_cum.48" u 1:(\$2) w l lt 2 lw 6   lc rgb "#377eb8" ti '\$\bm{^{48}}\$Ca p',\
"" u 1:(\$3) w l lt 3 lw 6  lc rgb "#377eb8" noti,\
"" u 1:(\$4) w l lt 4 lw 6  lc rgb "#377eb8" noti,\
"kin_cum.48" u 1:(\$5) w l lt 2 lw 6   lc rgb "#4daf4a" ti '\$\bm{^{48}}\$Ca n',\
"" u 1:(\$6) w l lt 3 lw 6  lc rgb "#4daf4a" noti,\
"" u 1:(\$7) w l lt 4 lw 6  lc rgb "#4daf4a" noti,\
"kin_cum.56" u 1:(\$2) w l lt 2 lw 6   lc rgb "#984ea3" ti '\$\bm{^{56}}\$Fe p',\
"" u 1:(\$3) w l lt 3 lw 6  lc rgb "#984ea3" noti,\
"" u 1:(\$4) w l lt 4 lw 6  lc rgb "#984ea3" noti,\
"kin_cum.56" u 1:(\$5) w l lt 2 lw 6   lc rgb "#ff7f00" ti '\$\bm{^{56}}\$Fe n',\
"" u 1:(\$6) w l lt 3 lw 6  lc rgb "#ff7f00" noti,\
"" u 1:(\$7) w l lt 4 lw 6  lc rgb "#ff7f00" noti








??

pdflatex $file
rm $file-inc-eps-converted-to.pdf
#rm $file-inc.eps
rm $file.{aux,log}
