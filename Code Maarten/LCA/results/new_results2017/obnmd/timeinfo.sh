#!/bin/bash

##
# this script goes through the files
# dens_ob4...
# in the folder you supply via argument
# and extracts the elapsed time (# elapsed time is: [....] s)
# and the mass number A (# A = [...] ...)
#
# usage:
#
# ./timeinfo.sh [inputfolder]
#
##


for file in $1/dens_ob4*;
do
    elapsedtime=$(cat ${file} | sed -n  's/# elapsed time is:\s*\(\S*\)\s*s.*/\1/p')
    A=$(cat ${file} | sed -n 's/#\s*A\s*=\s*\(\S*\).*/\1/p')
    Z=$(cat ${file} | sed -n 's/.*Z\s*=\s*\(\S*\).*/\1/p')
    printf "%4d  %4d   %8.2f\n" "${A}" "${Z}" "${elapsedtime}"
done | sort -n
