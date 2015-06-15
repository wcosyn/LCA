#!/bin/bash
#for file in `ls wscommomgsl_*.Pb`
#do
#  if [ `wc $file -l | cut -d' ' -f1` -ne 62 ]
#  then
#    echo `wc $file -l | cut -d" " -f1`
#  fi
#done

for file in `ls wscommomgsl_*.Pb.~*`
do
  base=${file%.\~?\~}

  if [ `ls $base* | wc -l` -eq 1 ]
  then
    mv -vi $file $base

  elif [ `ls $base* | wc -l` -eq 2 ]
  then
    file1=`ls $base* | head -n1`
    file2=`ls $base* | tail -n1`
    diff=`diff -q $file1 $file2`
    if [ -z "$diff" ]
    then
      rm -v $file2
    else
      ./integral.py $file1
      ./integral.py $file2
    fi

#    val1=`head -n2 $file1 | tail -n1 | cut -f2`
#    val2=`head -n2 $file2 | tail -n1 | cut -f2`
#
#    if [ $val1 == 0 -a $val2 != 0 ]
#    then
#      rm -v $file1
#    fi
#
#    if [ $val2 == 0 -a $val1 != 0 ]
#    then
#      rm -v $file2
#    fi
    
  else
    echo `ls $base*`
  fi
done

if [ -e integral.Pb ]
then
  rm integral.Pb
fi

for file in `ls wscommomgsl_*.Pb`
do
  ./integral.py $file >> integral.Pb
done


