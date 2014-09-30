#!/bin/bash
file="test.gp"
echo -e "" > $file

while read line 
do
  label=`echo $line | awk '{ print $1 }'`
  if [ $label != "#" ]
  then
    X=`echo $line | awk '{ print $2 }'`
    x=`echo $line | awk '{ print $4 }'`
    y=`echo $line | awk '{ print $3 }'`
    echo "set label '^{$label}$X' at $x,$y-0.04 " >> $file
  fi
done < EMC
