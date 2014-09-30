#!/usr/bin/bash

for file in `ls dens_ob4.*.111.*`
do
  echo $file
  (head -n 8 $file && head -n -3 $file | tail -n +9 | sort -n && tail -n 3 $file) > sorted_$file
  cp sorted_$file $file
done

#outfile=dens_ob.D
#array=(sorted_dens_ob.00.110.D.*)
#echo ${array[@]}
#echo "./merge.py $outfile ${array[@]}"
#./merge.py $outfile ${array[@]}

