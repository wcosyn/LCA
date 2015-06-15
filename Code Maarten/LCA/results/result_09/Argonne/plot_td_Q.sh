#!/usr/bin/bash


#python calc_norm_dens_td.py 0 000 dens_td.-11.He 1
#python calc_norm_dens_td.py 0 000 dens_td.11.He 3
#echo $NORM $NORM2 > norms_calc.txt
#./plot_td_Q.plt 000 he4 He 0 0 $NORM $NORM2

echo "" > norms_calc.txt

#ARGS x[fm] xxx[fm*100] index 
function calc
{
  NORM=`python calc_norm_dens_td.py $1 $2 dens_td.-11.111.He 1` 
  NORM2=`python calc_norm_dens_td.py $1 $2 dens_td.11.111.He 3`
  echo -e "$1 \t $NORM \t $NORM2 " >> norms_calc.txt
  ./plot_td_Q.plt $2 he4 He $1 $3 $NORM $NORM2
}

calc 0 000 0
calc 0.05 005 1
calc 0.10 010 2
calc 0.15 015 3
calc 0.20 020 4
calc 0.25 025 5
calc 0.30 030 6
calc 0.35 035 7
calc 0.40 040 8
calc 0.45 045 9
calc 0.50 050 10
calc 0.55 055 11
calc 0.60 060 12
calc 0.65 065 13
calc 0.70 070 14
calc 0.75 075 15
calc 1.00 100 20
calc 1.25 125 25
calc 1.50 150 30
calc 1.75 175 35
calc 2.00 200 40


#NORM=`python fit_norm_dens_td.py 0.05 005 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.05 005 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 005 he4 He 5 1 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.1 010 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.1 010 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 010 he4 He 10 2 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.15 015 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.15 015 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 015 he4 He 15 3 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.20 020 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.20 020 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 020 he4 He 20 4 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.25 025 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.25 025 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 025 he4 He 25 5 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.30 030 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.30 030 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 030 he4 He 30 6 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.35 035 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.35 035 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 035 he4 He 35 7 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.40 040 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.40 040 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 040 he4 He 40 8 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.45 045 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.45 045 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 045 he4 He 45 9 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.50 050 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.50 050 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 050 he4 He 50 10 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.55 055 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.55 055 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 055 he4 He 55 11 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.60 060 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.60 060 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 060 he4 He 60 12 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.65 065 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.65 065 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 065 he4 He 65 13 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.70 070 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.70 070 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 070 he4 He 70 14 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 0.75 075 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 0.75 075 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 075 he4 He 75 15 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 1 100 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 1 100 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 100 he4 He 100 20 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 1.25 125 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 1.25 125 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 125 he4 He 125 25 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 1.50 150 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 1.50 150 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 150 he4 He 150 30 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 1.75 175 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 1.75 175 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 175 he4 He 175 35 $NORM $NORM2
#NORM=`python fit_norm_dens_td.py 2 200 dens_td.-11.He 1`
#NORM2=`python fit_norm_dens_td.py 2 200 dens_td.11.He 3`
#echo $NORM $NORM2 >> norms.txt
#./plot_td_Q.plt 200 he4 He 200 40 $NORM $NORM2
