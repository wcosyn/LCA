#!/bin/bash

echo "DO NOT USE"
exit

if [ $# -ne 1 ]
then
	echo " One argument needed "
	exit
fi
newdir=$1
inputdir=`pwd`
cd $newdir
tmpfile="$inputdir/tmp.dat"
echo $tmpfile
for file in `ls recmosh*.dat`
do
	if [ -e $inputdir/$file ]
	then
		newfile=$file.new
		echo $newfile
		number=${file#recmosh}
		number=${number%.dat}
		lines=`wc $file | awk '{ print $1 }'`
		linesi=`wc $inputdir/$file | awk '{ print $1 }'`
		tail $file -n `expr $lines - 1` > $tmpfile  
		tail $inputdir/$file -n `expr $linesi - 1` >> $tmpfile  
		echo "$number" > $newfile
		cat $tmpfile | sort | uniq  >> $newfile
	else
		cp -vi $file $inputdir
	fi
done
cp -v *.dat.new $inputdir
cd $inputdir
for file in `ls *.dat.new`
do
	oldfile=${file%.new}
	mv $oldfile $oldfile.old
	mv -vi $file $oldfile
done





