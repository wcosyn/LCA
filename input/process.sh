#!/bin/bash


for file in `ls *.dat`
do
  if [  -e $file.\~1\~ ]
  then
    echo $file 
    if [ ! -e $file.\~2\~ ]
    then
      mv $file $file.\~2\~
    elif [ ! -e $file.\~3\~ ]
    then
      mv $file $file.\~3\~
    elif [ ! -e $file.\~4\~ ]
    then
      mv $file $file.\~4\~
    elif [ ! -e $file.\~5\~ ]
    then
      mv $file $file.\~5\~
    elif [ ! -e $file.\~6\~ ]
    then
      mv $file $file.\~6\~
    elif [ ! -e $file.\~7\~ ]
    then
      mv $file $file.\~7\~
    elif [ ! -e $file.\~8\~ ]
    then
      mv $file $file.\~8\~
    else
      echo " $file not renamed "
      exit
    fi 
  fi
done
echo "rename.. [OK]"

for file1 in `ls *.\~1\~`
do
	file=${file1%.\~1\~}
	if [ -e $file ]
	then
	  echo "ERR: $file and $file1 exists !!!"
	  exit
	fi
	tmpfile="tmp.dat"
#	echo $testfile
	number=${file#recmosh}
	number=${number%.dat}
#	echo $number
	file1=$file.\~1\~
	if [ -e $file1 ]
	then
		lines=`wc $file1 | awk '{ print $1 }'`
		tail $file1 -n `expr $lines - 1` > $tmpfile  
		mv $file1 trash/
	
	fi
	file2=$file.\~2\~
	if [ -e $file2 ]
	then
		lines=`wc $file2 | awk '{ print $1 }'`
		tail $file2 -n `expr $lines - 1` >> $tmpfile  
		mv $file2 trash/
	
	fi
	file3=$file.\~3\~
	if [ -e $file3 ]
	then
		lines=`wc $file3 | awk '{ print $1 }'`
		tail $file3 -n `expr $lines - 1` >> $tmpfile  
		mv $file3 trash/
	
	fi
	file4=$file.\~4\~
	if [ -e $file4 ]
	then
		lines=`wc $file4 | awk '{ print $1 }'`
		tail $file4 -n `expr $lines - 1` >> $tmpfile  
		mv $file4 trash/
	
	fi
	file5=$file.\~5\~
	if [ -e $file5 ]
	then
		lines=`wc $file5 | awk '{ print $1 }'`
		tail $file5 -n `expr $lines - 1` >> $tmpfile  
		mv $file5 trash/
	
	fi
	file6=$file.\~6\~
	if [ -e $file6 ]
	then
		lines=`wc $file6 | awk '{ print $1 }'`
		tail $file6 -n `expr $lines - 1` >> $tmpfile  
		mv $file6 trash/
	
	fi
	file7=$file.\~7\~
	if [ -e $file7 ]
	then
		lines=`wc $file7 | awk '{ print $1 }'`
		tail $file7 -n `expr $lines - 1` >> $tmpfile  
		mv $file7 trash/
	
	fi
	file8=$file.\~8\~
	if [ -e $file8 ]
	then
		lines=`wc $file8 | awk '{ print $1 }'`
		tail $file8 -n `expr $lines - 1` >> $tmpfile  
		mv $file8 trash/
	
	fi
	echo "$number" > $file
	cat $tmpfile | sort | uniq >> $file
done
		
	

#for file1 in `ls | grep \~1\~`
#do
#	echo $file1
#	file=${file1%.\~1\~}
#	file2=$file.\~2\~
#	if [ ! -e $file2 ]
#	then
#		echo "mv $file $file2"
#	fi
#done
