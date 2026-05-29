#!/bin/bash

if [ "$1" == "" ]
then
    echo "Min count (for both) needed"
    exit 1
fi
min_count=$1

for thresh in `cat fracs`
do
    echo $thresh
    rc="rc/gfga-${min_count}-${thresh}.rc"
    rc2="rc/tax_gfga-${min_count}-${thresh}.rc"
    gawk -v min_count=$min_count \
         'NR > 1 && $2 > min_count {print $1, $3/$2}' gf_vs_ga.dat | \
        gawk -v thresh=$thresh '$2 >= thresh {print $1, 3000}' > $rc
    echo "  $rc"
    rcmerge.sh -r -v2 tax.rc $rc > $rc2
    echo "  $rc2"
done
