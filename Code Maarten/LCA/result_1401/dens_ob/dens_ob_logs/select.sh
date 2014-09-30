#!/usr/bin/sh

for file in `ls dens_ob_*`
do
  echo $file
  egrep ^mf $file -A1 -B1 | tail -n2 
  echo 
done

