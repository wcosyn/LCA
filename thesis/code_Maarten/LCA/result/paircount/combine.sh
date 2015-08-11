#!/bin/bash

newfile="PairsNP"
echo -e "# NP Pairs \n# Norm: N*Z\n# A \t Z \t norm \t (n l S) = (- 0 -) \t (- 0 0) \t (- 0 1)\t (0 0 1) \t (n l S L) = (0 0 1 0) \t (- 0 1 0)" > $newfile
for file in `ls Pairs.*.-1.1`
do
  echo $file
  A=`head -n1 $file | cut -d " " -f2`
  Z=`head -n1 $file | cut -d " " -f3`
  norm=`grep "(- - -)"  $file | cut -d/ -f2`
  first=`grep '(- 0 -)' $file | cut  -f2`
  second=`grep '(- 0 0)' $file | cut  -f2`
  third=`grep '(- 0 1)' $file | cut  -f2`
  fourth=`grep '(0 0 1)' $file | cut  -f2`
  fourthb=`grep '(0 0 0)' $file | cut  -f2`
  fifth=`grep '(0 0 1 0)' $file | cut -f2`
  sixth=`grep '(- 0 1 0)' $file | cut -f2`
  echo -e "$A \t $Z \t $norm \t $first \t $second \t $third \t $fourth \t $fourthb \t $fifth \t $sixth" >> $newfile
done


newfile="PairsPP"
echo -e "# PP Pairs \n# Norm: Z(Z-1)/2\n# A \t Z \t norm \t (n l S) = (- 0 -) \t (- 0 0) \t (0 0 0) \t (n l S L) = (0 0 0 0)" > $newfile
for file in `ls Pairs.*.1.1`
do
  A=`head -n1 $file | cut -d " " -f2`
  Z=`head -n1 $file | cut -d " " -f3`
  norm=`grep "(- - -)"  $file | cut -d/ -f2`
  first=`grep '(- 0 -)' $file | cut  -f2`
  second=`grep '(- 0 0)' $file | cut  -f2`
  third=`grep '(0 0 0)' $file | cut  -f2`
  fourth=`grep '(0 0 0 0)' $file | cut -f2`
  echo -e "$A \t $Z \t $norm \t\t $first \t $second \t $third \t\t $fourth" >> $newfile
done

newfile="PairsNN"
echo -e "# NN Pairs \n# Norm: N(N-1)/2\n# A \t Z \t norm \t (n l S) = (- 0 -) \t (- 0 0) \t (0 0 0) \t (n l S L) = (0 0 0 0)" > $newfile
for file in `ls Pairs.*.-1.-1`
do
  A=`head -n1 $file | cut -d " " -f2`
  Z=`head -n1 $file | cut -d " " -f3`
  norm=`grep "(- - -)"  $file | cut -d/ -f2`
  first=`grep '(- 0 -)' $file | cut  -f2`
  second=`grep '(- 0 0)' $file | cut  -f2`
  third=`grep '(0 0 0)' $file | cut  -f2`
  fourth=`grep '(0 0 0 0)' $file | cut -f2`
  echo -e "$A \t $Z \t $norm \t\t $first \t $second \t $third \t\t $fourth" >> $newfile
done
