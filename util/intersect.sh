#!/bin/bash

if [ "$2" == "" ]
then
    echo "intersect.sh <react-cost-fn> <scenario-fn>"
    echo "    Extracts scenario from <scenario-fn> where"
    echo "    reactions are in <react-cost-fn>"
    echo "    Replaces scenario cost with react-cost cost"
    exit 1
fi

reacts=$1
dat=$2

echo "# Intersect scenario $data with reactions $reacts"
echo "# On `date`"
grep '^MET' $dat
# Get biomass react and cost-1 reactions
gawk '/^REACT/ && $3 <= 1' $dat
grep ^REACT $dat | kjoin $reacts - -j2 2 | gawk '{$5 = $2; $1 = ""; $2 = ""; print;}' | sed 's/^\s*//'
