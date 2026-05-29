#!/bin/bash

# - Clean up the names of reactions,
# - add reverse reacts
# - remove biomass reactions
sed -e 's/"//g' -e 's/^R_//' $1 | gawk '$1 !~ /[Bb]iomass/ {print $1, $2 "\n" $1 "__rev", $2}'
