#!/bin/bash
newfile=wscommomgsl.Al
tmpfile=tmp.Al
if [ -e $newfile ]
then
	mv $newfile $newfile.old -v --backup=t
fi



for file in `ls wscommomgsl_*.Al`
do
  shell1=${file:12:3}
  tmp=${file##*_}
  shell2=${tmp:0:3}
  if [ "$shell1" = "025" -o "$shell2" = "025" ]
  then
    if [ "$shell1" = "$shell2" ]
    then
      factor=0.66666667
    else
      factor=0.83333333
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


