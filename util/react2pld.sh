#!/bin/bash

Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo "Usage: `basename $0` react_fn univ.pld [out.pld]" >& 2
    echo "       Extracts METs and REACTs from a PLD file for a" >&2
    echo "       given set of reactions." >&2
    echo "  Switches" >&2
    echo "    " >&2
    echo "  Args" >&2
    echo "     react_fn: list of reaction names " >&2
    echo "     univ.pld: Plum pld file with all METs and REACTs" >&2
    echo "      out.pld: output filename [out.pld]" >&2
    exit 1
}

a=0
while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-a")  a="$1" ; shift ;;
      "-?")  Usage ;;
  "--help")  Usage ;;
      *) Usage "Unknown option" $arg ;;
  esac
done

# Usage 
if [ $# -lt 2 ]
then
    Usage
fi
react_fn="$1"
univ_fn="$2"
out_fn="${3:-out.pld}"

if [ ! -f $react_fn ] 
then
    Usage "Couldn't open file $react_fn" 
fi
if [ ! -f $univ_fn ] 
then
    Usage "Couldn't open file $univ_fn" 
fi

tmp="/tmp/react2pld.tmp"
tmp2="/tmp/react2pld2.tmp"

grep ^REACT $univ_fn | kjoin $react_fn - -j2 2 -o2 > $tmp
gawk '/^REACT/ {for (k = 5; k <= NF; k += 2) {print $k;}}' $tmp | sort | uniq > $tmp2

echo "# Created by react2pld.sh from $react_fn and $univ_fn on `date`" > $out_fn
# Grab MET lines
grep ^MET $univ_fn | kjoin $tmp2 - -j2 2 -o2 >> $out_fn
# Grab REACT lines
cat $tmp >> $out_fn

met_count=`cat $tmp2 | wc -l`
react_count=`cat $react_fn | wc -l`
out_count=`cat $tmp | wc -l`
echo "Done. Output in $out_fn"
echo "  $met_count metabolites and $out_count reactions"
if [ $react_count != $out_count ]
then
    echo "NOTE: $react_count lines in $react_fn, but only $out_count matches found"
    grep ^REACT $univ_fn | kjoin $react_fn - -j2 2 -r1 > $tmp
    echo "View missing reactions in $tmp"
else
   echo "All $react_count lines in $react_fn found a match"
fi       
