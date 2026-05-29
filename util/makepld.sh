#!/bin/bash
# Run a command
function run() {
    echo "$@"
    "$@" 
}
# Run a command with redirect to first arg
function run_redirect() {
    out_file=$1
    shift
    echo "$@ > $out_file"
    "$@" > $out_file
    if [ $? -ne 0 ]
    then
        pause
    fi
}
# Display message and wait for a button-press
function pause() {
 read -s -n 1 -p "Press any key to continue . . ."
 echo ""
}

out=out.pld

if [ "$2" == "" ]
then
    echo "Usage: makepld <all.pld> <flux_file> [out.pld]"
    echo "    all.pld: Our MET/REACT file with all reactions"
    echo "  flux_file: list of <reaction> "
    echo "    "
    echo "           out.pld: output [$out]"
    exit 1
fi

all_pld=$1
flux=$2
if [ "$3" != "" ]
then
    out=$3
fi
tmp_react="/tmp/makepld1.tmp"
tmp_met="/tmp/makepld2.tmp"

tmp3="/tmp/makepld3.tmp"
tmp4="/tmp/makepld4.tmp"

# Select reactions
echo "Select reactions with known flux"
getcol 1 $flux | uniq | kjoin -c $all_pld - -j1 2 -o1 > $tmp_react

echo "Check we know them all"
getcol 1 $flux | uniq | kjoin -c $all_pld - -j1 2 -o1 -r2 > $tmp3
if [ `wc -l < $tmp3` -gt 0 ]
then
    echo "Some reactions not known:"
    cat $tmp3
fi

# Which compounds are used?
# Select used metabolites
echo "Select metabolites from known-flux reactions"
gawk '/^REACTION/ {for (i=5;i<NF;i+=2) {if ($i == "#") {break;} print $i;}}' $tmp_react | sort | uniq > $tmp_met


echo "# Created by makepld `date`" > $out
echo "# Using makepld $flux $all_pld" >> $out

# Grab mets
echo "Write metabolites"
grep '^MET' $all_pld | kjoin $tmp_met - -j2 2 -o2 >> $out 

cat $tmp_react >> $out

echo "Done. Output in $out"



