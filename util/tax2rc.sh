#!/bin/bash


tax="$1"
col="$2"
out=${3:-"tax_raw.rc"}

if [ "$2" == "" ]
then
    echo '
Usage tax2rc.sh <taxonomy-csv> <col-num> [<out-fn>]"
    Takes csv containing the taxonomy weights for all organisms, and
    extracts the weights (reaction costs) for a single one.
    - removes quotes from reaction name
    - removes "R_" from start of reaction name
    - excludes biomass reactions
    - duplicates for __rev

Args
    taxonomy-csv: csv filename of taxonomy weights [none]
         col-num: the 1-based col num of weight col to extract [none]
          out-fn: filename for output [$out]
'
    exit 1
fi

echo "# Extracted from $tax by tax2rc.sh on `date`" > $out
name=`head -1 "$tax" | getcol -F, $col -`
echo "Exracting weights for $name"
getcol -F, 1 $col "$tax" | sed -e 's/"//g' -e 's/^R_//' | grep -v biomass | gawk '
NR > 1 {print $1, $2; print $1 "__rev", $2}' >> $out
count=`wc  -l < $out`

echo "Done. $count weights extracted. Output in $out"
