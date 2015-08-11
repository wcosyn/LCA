#!/usr/bin/bash
# no ARG

file="dens_rel_ct"

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

set term epslatex standalone color  ",bx," 17 dl 3  size xsize,ysize header '\usepackage{bm}'
set output "$file".".tex"

hbarc = 0.197327

#set style line 1 lt 1 lw 6 lc rgb "#000000"
set style line 2 lt 1 lw 6 lc rgb "#E41A1C"
set style line 3 lt 2 lw 6 lc rgb "#377eb8"
set style line 4 lt 5 lw 6 lc rgb "#4daf4a"
set style line 5 lt 4 lw 6 lc rgb "#984EA3" 
set style line 6 lt 7 lw 6 lc rgb "#ff7f00"
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
set xr[0:4]
set xtics 0,1,5 offset -0.3
set mytics 5
set yr[1e-4:2e1]
set ytics 1e-4, 10, 1e3 offset 0, -0.25
set logscale y
set key samplen 1.9
set key top right height 1 Right spacing 1.05

set pointsize 1.2


#start plotting
set multiplot

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(2)
set bmargin at screen bot(2)


set ylabel '\$\bm{n^{[2]}(k_{12})}\$ [fm\$\bm{^3}\$]' offset 2.1, 0
set format y "\$\\\bm{10^{%L}}\$"
set xlabel ''
set format x ''
#set xlabel "\$\\\bm{k_{12}}\$ [fm\$\\\bm{^{-1}}\$]"
#set format x "%1.0f"

A=4
X="He"

plot  "./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$5)/A/(A-1)*2 w l ls 2 ti '\$\bm{^{'.A.'}}\$'.X.' LCA',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$2)/A/(A-1)*2 w l ls 3 ti 'IPM',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$6)/A/(A-1)*2 w l ls 4 ti 'cen',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$7)/A/(A-1)*2 w l ls 5 ti 'ten'

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(2)
set bmargin at screen bot(2)

#set ylabel '\$\bm{n^{[2]}(k_{12})}\$ [fm\$\bm{^3}\$]' offset 2.1, 0
#set format y "\$\\\bm{10^{%L}}\$"
set ylabel ''
set format y ''
set xlabel ''
set format x ''
#set xlabel "\$\\\bm{k_{12}}\$ [fm\$\\\bm{^{-1}}\$]"
#set format x "%1.0f"

A=9
X="Be"


plot  "./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$5)/A/(A-1)*2 w l ls 2 ti '\$\bm{^{'.A.'}}\$'.X.' LCA',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$2)/A/(A-1)*2 w l ls 3 ti 'IPM',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$6)/A/(A-1)*2 w l ls 4 ti 'cen',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$7)/A/(A-1)*2 w l ls 5 ti 'ten'

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(3)
set rmargin at screen right(3)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(2)
set bmargin at screen bot(2)

#set ylabel '\$\bm{n^{[2]}(k_{12})}\$ [fm\$\bm{^3}\$]' offset 2.1, 0
#set format y "\$\\\bm{10^{%L}}\$"
set ylabel ''
set format y ''
set xlabel ''
set format x ''
#set xlabel "\$\\\bm{k_{12}}\$ [fm\$\\\bm{^{-1}}\$]"
#set format x "%1.0f"

A=12
X="C"


plot  "./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$5)/A/(A-1)*2 w l ls 2 ti '\$\bm{^{'.A.'}}\$'.X.' LCA',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$2)/A/(A-1)*2 w l ls 3 ti 'IPM',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$6)/A/(A-1)*2 w l ls 4 ti 'cen',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$7)/A/(A-1)*2 w l ls 5 ti 'ten'

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set ylabel '\$\bm{n^{[2]}(k_{12})}\$ [fm\$\bm{^3}\$]' offset 2.1, 0
set format y "\$\\\bm{10^{%L}}\$"
#set ylabel ''
#set format y ''
#set xlabel ''
#set format x ''
set xlabel "\$\\\bm{k_{12}}\$ [fm\$\\\bm{^{-1}}\$]"
set format x "%1.0f"

A=16
X="O"


plot  "./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$5)/A/(A-1)*2 w l ls 2 ti '\$\bm{^{'.A.'}}\$'.X.' LCA',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$2)/A/(A-1)*2 w l ls 3 ti 'IPM',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$6)/A/(A-1)*2 w l ls 4 ti 'cen',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$7)/A/(A-1)*2 w l ls 5 ti 'ten'

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(1)
set bmargin at screen bot(1)

#set ylabel '\$\bm{n^{[2]}(k_{12})}\$ [fm\$\bm{^3}\$]' offset 2.1, 0
#set format y "\$\\\bm{10^{%L}}\$"
set ylabel ''
set format y ''
#set xlabel ''
#set format x ''
set xlabel "\$\\\bm{k_{12}}\$ [fm\$\\\bm{^{-1}}\$]"
set format x "%1.0f"

A=27
X="Al"


plot  "./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$5)/A/(A-1)*2 w l ls 2 ti '\$\bm{^{'.A.'}}\$'.X.' LCA',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$2)/A/(A-1)*2 w l ls 3 ti 'IPM',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$6)/A/(A-1)*2 w l ls 4 ti 'cen',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$7)/A/(A-1)*2 w l ls 5 ti 'ten'

#--------------------------------------------------
#  set horizontal margins for first column
set lmargin at screen left(3)
set rmargin at screen right(3)
#  set horizontal margins for lasts row (top)
set tmargin at screen top(1)
set bmargin at screen bot(1)

#set ylabel '\$\bm{n^{[2]}(k_{12})}\$ [fm\$\bm{^3}\$]' offset 2.1, 0
#set format y "\$\\\bm{10^{%L}}\$"
set ylabel ''
set format y ''
#set xlabel ''
#set format x ''
set xlabel "\$\\\bm{k_{12}}\$ [fm\$\\\bm{^{-1}}\$]"
set format x "%1.0f"

A=40
X="Ca40"
Xkey="Ca"


plot  "./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$5)/A/(A-1)*2 w l ls 2 ti '\$\bm{^{'.A.'}}\$'.Xkey.' LCA',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$2)/A/(A-1)*2 w l ls 3 ti 'IPM',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$6)/A/(A-1)*2 w l ls 4 ti 'cen',\
"./dens_rel.00.110.".X.".-1-1-1-1.-1-1" u 1:(\$7)/A/(A-1)*2 w l ls 5 ti 'ten'

unset multiplot



??

pdflatex $file
rm $file-inc-eps-converted-to.pdf
#rm $file-inc.eps
rm $file.{aux,log}
