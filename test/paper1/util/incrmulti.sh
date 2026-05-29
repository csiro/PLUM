#!/bin/bash
function run() {
    echo "$@"
    eval "$@"
}
function die()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    exit 1
}


dir="multi"
if [ ! -d $dir ]
then
    mkdir -v $dir
fi

if [ "$2" == "" ]
then
    die "pld filename and sdfn-lis required"
fi
pld="$1"
sd_lis="$2"

res="incrmulti.res"
tmp=/tmp/incrmulti.tmp

allsols="${dir}/sol.lis"

plum_args="-v incr -rc tax.rc"

echo -n "" > $allsols

for sd in `cat $sd_lis`
do
    b=`basename $sd .sd`
    sol="${dir}/${b}.sol"
    echo $sol >> $allsols

    if [ ! -f $sol ]
    then
        run "plum $plum_args -sd $sd -o $sol -O $res $pld -C $cost"
    else
        echo "Already have sol for $sd"
    fi
done

reacts=$dir/reacts.lis

# Also grab gene-indicated reactions
gawk '/^REACT/ && $3 == 1 {print $2, $3}' $pld > $tmp

cat $tmp `cat $allsols` | getcol -c 1 | sort | uniq > $reacts

model=incrmulti.pld
react2pld.sh $reacts $pld $model

echo "Done. Model in $model"


