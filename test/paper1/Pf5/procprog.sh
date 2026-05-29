#!/bin/bash

if [ "$1" == "" ]
then
    echo "Progress fn required"
    exit 1
fi

prog=$1
tmp="/tmp/procprog.tmp"

nc=7
grep Method $prog > $tmp
getcol iter weight+1 weight+2 weight+3 weight+4 weight+5 weight+6 weight+7 $prog >> $tmp


splitcols -N $tmp | cat -s 

echo "" 
echo "L Sum"
getcol iter weight+1 weight+2 weight+3 weight+4 weight+5 weight+6 $prog | gawk '
/^[0-9]/ {sum=0; for (k=2; k<= NF; k++) {sum += $k;} print $1, sum}' 

