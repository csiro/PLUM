#!/bin/bash
mult=1.0
if [ $# -ge 1 ]
then
    mult=$1
    shift
fi
sd_list="growth.lis nogrowth.lis"
if [ $# -ge 1 ]
then
    sd_list="$@"
fi
#echo "sd list is $sd_list"
tmp=/tmp/biolog.tmp
echo -n "" > $tmp
for sd in `cat $sd_list`
do
	grep -H biolog $sd | sed 's/:biolog//' >> $tmp
done 

sort -nr -k 2,2 $tmp | gawk -v "mult=$mult" 'NR == 1 {mx=$2} {print $1, mult*$2/mx}'
