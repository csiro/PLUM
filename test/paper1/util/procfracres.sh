#!/bin/bash

if [ "$1" == "" ]
then
    echo "Fracs fn required"
    exit 1
fi
tmp="/tmp/procfracres.tmp"
grep biolog sd/*.sd | sed 's/:biolog//' > $tmp
for f in `cat "$1"`
do
    echo "--------------------------"
    echo -n $f
    dir="testall-$f"
    bio="$dir/biomass.$f"
    getcol sd_fn biomass_flux $dir/testall.res > $bio
    echo -n "  Num no biomass "
    gawk '$2 == 0' $bio | wc -l
    echo -n "Cost of union sol: "
    getcol -c 3 $dir/union.txt | kstat -q -show sum
    echo "Different to base..."
    echo "Base $f biolog"
    gawk '$2 == 0' $bio | kjoin biomass.base - | kjoin - $tmp | getcol 1 2 4 6 | gawk '$2 != $3'
done

