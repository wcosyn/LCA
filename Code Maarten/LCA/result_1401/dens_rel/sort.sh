#!/usr/bin/bash

for file in `ls dens_rel.00.110.$1.*`
do
  echo $file
  (head -n 8 $file && head -n -4 $file | tail -n +9 | sort -n && tail -n 4 $file) > sorted_$file
done

#outfile=dens_ob.D
#array=(sorted_dens_ob.00.110.D.*)
#echo ${array[@]}
#echo "./merge.py $outfile ${array[@]}"
#./merge.py $outfile ${array[@]}

