#!/bin/bash

if [ "$1" == "" ]
then
    echo "Fracs fn required"
    exit 1
fi
for f in `cat "$1"`
do
    echo $f
    rc="rc/tax_gfga-${f}.rc"
    if [ ! -f $rc ]
    then
        echo "React cost file not found: $rc"
        exit 1
    fi
    testall.sh -a "react_cost_fn=$rc max_react_cost=2001" Pf5_agora2.pld sd.lis
    if [ -d testall-${f} ]
    then
        rm -r testall-${f}
    fi
    mv testall testall-${f}
done
