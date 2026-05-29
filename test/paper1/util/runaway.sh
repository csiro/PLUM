#!/bin/bash
# Run a command
function run() {
    echo "$@"
    eval "$@"
}
function die()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    exit 1
}

dir="working"
if [ ! -d $dir ]
then
    mkdir -v $dir
fi

count=0
pld=$dir/runaway${count}.pld
if [ "$1" == "-c" ]
then
    shift
    count="$1"
    shift
    pld=$dir/runaway${count}.pld
    if [ ! -f $pld ]
    then
        die "Lost file $pld"
    fi
elif [ "$1" == "" ]
then
    die "pld filename required"
else
    if [ ! -f $1 ]
    then
        die "Lost file $1"
    fi
    cp -v $1 $pld
fi

plum_args="-q -sd sd/sucr.sd -rc tax.rc"
out="runaway.out"
sol="runaway.sol"

if [ $count -eq 0 ]
then
    echo -n "" > $out
fi

while true
do
    if [ -f $sol ]; then rm $sol; fi

    run plum $plum_args -o $sol $pld

    if [ ! -f $sol ]
    then
        die "No solution produced .. give up"
    fi
    
    top=( `sort -nr -k 2,2 -k 3,3 $sol | getcol -c 1 2 | head -1` )
    name=${top[0]}
    flux=${top[1]}

    if [ "$name" == "" ]
    then
        echo "No flux - abort!"
        break
    fi

    clear

    echo "Most active reaction is $name with flux $flux"

    iflux=`kcalc $flux rounddn`
    if [ $iflux -lt 10 ]
    then
        echo "Max flux is now $flux - finished!"
        break
    fi
    
    echo "Remove reaction $name"
    echo "$name $flux">> $out

    count=$(( count + 1 ))
    new_pld=$dir/runaway${count}.pld
    grep -v "^REACTION $name" $pld > $new_pld
    pld=$new_pld

done

echo "Done. `wc -l < $out` reactions removed."
echo "Last pld in $pld"
echo "List of removed reactions in $out"
