#!/bin/bash


pld=$1
sd_dir=${2:-"sd"}
sol_dir=${3:-"sol"}

if [ "$1" == "" ]
then
    echo '
Usage doall.sh <pld> [<sd_dir> [<sol-dir>]]
    Runs plum on given pld file for all carbons source (.sd) files
    in directory <sd_dir>. Put sols in <sol-dir>.
'
    exit 1
fi

if [ ! -d $sol_dir ]
then
    mkdir -v $sol_dir
fi

for sd in $sd_dir/*.sd
do
    echo "Running sd $sd"
    base=`basename $sd .sd`
    out="${sol_dir}/${base}.sol"
    # Remove old result
    grep -v "sd_fn $sd" res > res.tmp
    mv res.tmp res

    cmd="plum $pld -sd $sd -O res -o $out"
    echo "$cmd"
    $cmd
done
