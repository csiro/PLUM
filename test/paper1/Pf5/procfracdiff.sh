#!/bin/bash
function die() {
    echo "$@"
    exit 1
}

if [ "$1" == "" ]
then
    die "Fracs fn required"
fi

tmp="/tmp/procfracdiff.tmp"
grep biolog sd/*.sd | sed 's/:biolog//' > $tmp
tmp2="/tmp/procfracdiff2.tmp"

for frac in `cat "$1"`
do
    echo "=========================="
    echo -n $frac
    dir="testall-$frac"
    bio="$dir/biomass.$frac"
    getcol sd_fn biomass_flux $dir/testall.res > $bio
    
    echo -n "  Num no biomass "
    gawk '$2 == 0' $bio | wc -l
    echo -n "Cost of union sol: "
    getcol -c 3 $dir/union.txt | kstat -q -show sum
    echo "Different to base..."
    echo "Base $frac biolog"
    gawk '$2 == 0' $bio | kjoin biomass.base - | kjoin - $tmp | getcol 1 2 4 6 | sed -e 's|sd/||' -e 's|\.sd||' | gawk '$2 != $3 && $4 >= 0.3' | tee $tmp2

    for cs in `getcol 1 $tmp2`
    do
        echo "--------------------------"
        echo "Reacts in base for $cs forbidden in $frac"
        base="testall-base/$cs.sol"
        if [ ! -f $base ]; then die "Lost base sol for $cs: $base"; fi

        rc="rc/gfga-$frac.rc"
        if [ ! -f $rc ]; then die "Lost react cost for $frac: $rc"; fi

        echo "React Cost Flux"
        getcol -c 1 3 2 $base | gawk '$2 > 1' | kjoin - $rc -o1 | sort -k3,3n | kjoin - Pf5_agora2.pld -j2 2 
    done
    echo "--------------------------"
done

