#!/bin/bash
newfile=wscommomgsl.Pb
tmpfile=tmp.Pb
if [ -e $newfile ]
then
	mv $newfile $newfile.old -v --backup=t
fi

for file in `ls wscommomgsl_*.Pb`
do
	if [ ! -e $newfile ]
	then
		cp $file $newfile -v
	else	
		mv $newfile $tmpfile
		echo "./merge.py $file $tmpfile $newfile"
		./merge.py $file $tmpfile $newfile
	fi
done


