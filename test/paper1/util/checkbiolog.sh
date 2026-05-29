#!/bin/bash

Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo -n "Usage: `basename $0` [-r fn] res biolog-fn" >& 2
    echo '
  Check the biomass flux in a res file against biolog reference bands

  Switches
      -r: reference sd - used to calc translation mult [sd/sucr.sd]
  Args
       res: results file
       biolog-fn: reference biolog score file.
                  Format should be
                  <sd_fn> <min-val> <max_val> <ave-val>
' >&2
    exit 1
}

num_vals=3
ref_sd="sd/sucr.sd"
# Another opt method
while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-n")  num_vals="$1" ; shift ;;
      "-r")  ref_sd="$1" ; shift ;;
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

res="$1"
biolog="$2"

if [ ! -r $res ]
then
    Usage "Couldn't open file $res"
fi
if [ ! -r $biolog ]
then
    Usage "Couldn't open file $biolog"
fi

tmp1=/tmp/checkbiolog1.tmp
tmp2=/tmp/checkbiolog2.tmp
tmp3=/tmp/checkbiolog3.tmp
tmp4=/tmp/checkbiolog4.tmp
tmp5=/tmp/checkbiolog5.tmp

getcol sd_fn biomass_flux $res | kjoin - $biolog | sort -k 6,6n > $tmp1

grep $ref_sd $tmp1 > $tmp2
if [ $? != 0 ]
then
    echo "Referenece sd fn $ref_sd not found in $res / $biolog"
    exit 1
fi
ref_targ=`getcol 2 $tmp2`
ave_val=`getcol 6 $tmp2`
echo "Ref targ ($ref_sd) is $ref_targ, biolog ave $ave_val"
mult=`kcalc $ref_targ $ave_val /`
echo "Mult is $mult"

dig="checkbiolog.dig"
echo "H \"Biomass / Biolog comparison, using conversion mult $mult\"" > $dig
echo "T \"Created by checkbiolog.sh on `date`\"" >> $dig
echo "L \"Res file\"" >> $dig
gawk '{print "XC", $1, $2}' $tmp1 >> $dig
gawk '{printf ("LP %d %lf", NR-1, $2); for (k=7; k<=NF;k++){printf (" %s", $k);} print ""}' $tmp1 >> $dig
echo "" >> $dig
echo "L \"Biolog min\"" >> $dig
gawk -v mult=$mult '{print "XC", $1, $4 * mult}' $tmp1 >> $dig

echo "" >> $dig
echo "L \"Biolog max\"" >> $dig
gawk -v mult=$mult '{print "XC", $1, $5 * mult}' $tmp1 >> $dig

echo "" >> $dig
echo "L \"Ref Point\"" >> $dig
gawk -v mult=$mult '{print "XC", $1, $6 * mult}' $tmp2 >> $dig

echo "" >> $dig
echo "S 0" >> $dig
echo "L \"No Growth\"" >> $dig
left=`head -1 $tmp1| getcol 1`
right=`tail -1 $tmp1| getcol 1`
echo "XC $left 0.3" >> $dig
echo "XC $right 0.3" >> $dig

echo "Done Output in $dig. View with"
echo "dig $dig &"

