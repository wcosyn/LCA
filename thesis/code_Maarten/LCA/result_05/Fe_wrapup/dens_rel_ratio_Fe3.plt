set term epslatex standalone color  "default, 14" dl 3.5
#set term postscript enhanced color "Helvetica, 25"
#set term postscript enhanced color "Computer-Modern-Roman, 30"

set style line 2 lt 2 lw 6 lc rgb "#E41A1C"


set style line 3 lt 3 lw 6 lc rgb "#377eb8"
set style line 4 lt 4 lw 6 lc rgb "#4daf4a"
set style line 5 lt 5 lw 6 lc rgb "#984EA3"
set style line 6 lt 6 lw 6 lc rgb "#ff7f00"
set style line 7 lt 7 lw 6 lc rgb "#a65628"
set style line 8 lt 8 lw 6 lc rgb "#ffff33"
set style line 9 lt 9 lw 6 lc rgb "#999999"

set bmargin 3.2
set lmargin 6.8
set tmargin 0.3
set rmargin 1.2
set style data points
set border lw 3
set tics scale 1,0.5
#set key font ",25"

hbarc = 0.197327

#set xr[0:3.5*hbarc]
set xr[0:0.7]
set xlabel "$k_{12}$ [GeV]"
set yr[0:1.1]
set nokey


#set ylabel "n_2^{2n+l}(k) / n_2^{all}(k)" offset 1.5, 0
set ylabel "$n_2^{2n+l}(k_{12}) / n_2(k_{12})$" offset 1.5, 0
set key samplen 1.5
set key top left height 1



set object 1 rectangle behind from first 0,0 to first 0.3,1.1 fc rgb"light-gray" fs noborder


plot "dens_rel_all.E0.Fe3" u ($1)*hbarc:6 ls 2 w l ti "2n+l = 0" ,\
"dens_rel_all.E1.Fe3" u ($1)*hbarc:6 ls 3 w l ti "2n+l = 1" ,\
"dens_rel_all.E2.Fe3" u ($1)*hbarc:6 ls 4 w l ti "2n+l = 2" ,\
"dens_rel_all.E3.Fe3" u ($1)*hbarc:6 ls 5 w l ti "2n+l = 3" ,\
"dens_rel_all.E4.Fe3" u ($1)*hbarc:6 ls 6 w l ti "2n+l = 4" ,\
"dens_rel_all.E5.Fe3" u ($1)*hbarc:6 ls 7 w l ti "2n+l = 5" 

