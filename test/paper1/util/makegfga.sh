#!/bin/bash

rc_fn="$1"
list="gfga-10-1.0.rc gfga-10-0.9.rc gfga-10-0.8.rc" 

if [ "$1" == "" ] || [ ! -f $rc_fn ]
then
    if [ "$1" != "" ] 
    then
        echo "Base file $rc_fn not found"
    fi
    echo "Usage: makegfga.sh base.rc"
    echo "  Merges a base react cost file with a "
    echo "  list of gfga react cost files."
    echo "  The list is"
    echo "    $list"
    exit 1
fi

base="`echo $rc_fn | sed 's/\.rc$//'`"

for gfga in $list
do
    fn="../common/rc/$gfga"
    if [ ! -f $fn ]
    then
        echo "React cost file $fn not found"
        exit 1
    fi
    out="${base}-${gfga}"
    cmd="rcmerge.sh -r -v2 $rc_fn $fn > $out"
    echo $cmd
    eval $cmd
done
