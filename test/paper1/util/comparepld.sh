#!/bin/bash

pld1=$1
pld2=$2

if [ "$pld2" == "" ]
then
    echo "Two .pld file names required"
    exit 1
fi

tmp1="/tmp/comparepld1.tmp"
tmp2="/tmp/comparepld2.tmp"
tmp3="/tmp/comparepld3.tmp"
out1="/tmp/pld1only.tmp"
out2="/tmp/pld2only.tmp"

grep ^REACT $pld1 > $tmp1
grep ^REACT $pld2 > $tmp2

n1=`wc -l < $tmp1`
n2=`wc -l < $tmp2`
kjoin $tmp1 $tmp2 -j 2 > $tmp3
incommon=`wc -l < $tmp3`

kjoin $tmp1 $tmp2 -j 2 -r1 | sort -k3,3n > $out1
kjoin $tmp1 $tmp2 -j 2 -r2 | sort -k3,3n > $out2
n1only=`wc -l < $out1`
n2only=`wc -l < $out2`

pc1=`kcalc $incommon $n1 / 100 x round`
pc2=`kcalc $incommon $n2 / 100 x round`

pco1=`kcalc $n1only $incommon / 100 x round`
pco2=`kcalc $n2only $incommon / 100 x round`

printf "%4d %s\n" $n1 "Reactions in $pld1"
printf "%4d %s\n" $n2 "Reactions in $pld2"
printf "%4d %s\n" $incommon "Reactions in common (${pc1}%/${pc2}%)"
printf "%4d %s\n" $n1only "Reactions in $pld1 only (${pco1}%) - $out1"
printf "%4d %s\n" $n2only "Reactions in $pld2 only (${pco2}%) - $out2"

