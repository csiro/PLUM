#!/bin/bash

if [ "$2" == "" ]
then
    echo '
Usage mergepld.sh pld1 pld2 [out.pld]
   Merges two PLD files.
   - Keeps union of metabolites
   - Keeps lower cost reactions
'
    exit 1
fi

pld1="$1"
pld2="$2"
out="${3:-out.pld}"
tmp="/tmp/mergepld.tmp"

echo "# Created by mergepld.sh on `date`" > $out
echo "# by merging $pld1 and $pld2" >> $out
echo "# $pld1:" >> $out
grep "^#" < "$pld1" >> $out
echo "# $pld2:" >> $out
grep "^#" < "$pld2" >> $out

grep "^MET" "$pld1" > $tmp
cat $tmp >> $out

grep ^MET "$pld2" | kjoin - $tmp -j 2 -r1 >> $out

cat "$pld1" "$pld2" | grep "^REACT" | sort | choosemin -k 2 -n 1 -c 3 >> $out


echo "Done. Output in $out"
