#!/bin/bash
newfile=wscommomgsl.Fe
tmpfile=tmp.Fe
if [ -e $newfile ]
then
	mv $newfile $newfile.old -v --backup=t
fi

for file in `ls wscommomgsl_*.Fe`
do
  shell1=${file:12:3}
  tmp=${file##*_}
  shell2=${tmp:0:3}
  if [ "$shell1" = "037" -o "$shell2" = "037" ]
  then
    if [ "$shell1" = "$shell2" ]
    then
      factor=0.53571428571428571428
    else
      factor=0.75
    fi
  else
    factor=1
  fi

  if [ ! -e $newfile ]
  then
    cp $file $newfile -v
  else	
    mv $newfile $tmpfile
    echo "./merge.py $file $tmpfile $newfile $factor"
    ./merge.py $file $tmpfile $newfile $factor
  fi
  wait
done
