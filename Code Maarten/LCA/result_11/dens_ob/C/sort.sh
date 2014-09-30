#!/usr/bin/bash

for file in `ls dens_ob.00.110.C.*`
do
  echo $file
  (head -n 6 $file && head -n -3 $file | tail -n +7 | sort -n && tail -n 3 $file) > sorted_$file
done

outfile=dens_ob.C
array=(sorted_dens_ob.00.110.C.*)
echo ${array[@]}
echo "./merge.py $outfile ${array[@]}"
./merge.py $outfile ${array[@]}

