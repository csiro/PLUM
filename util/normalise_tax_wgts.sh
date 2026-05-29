#!/bin/bash

min_cost=1.5
max_cost=1999

Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo -n "Usage: `basename $0` [-r] rc-raw" >& 2
    echo '
       Normalise the reaction cost scores to values between ' $min_cost ' and ' $max_cost '
  Switches
     -r: Scores are likelihood, not cost (i.e. 0 in file is LEAST likely)

  Args
     rc-raw: file with <reaction-name> <raw-cost> pairs [tax_raw.rc]

' >&2
    exit 1
}

r_flag=0

while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-r")  r_flag=1 ;;
      "-?")  Usage ;;
  "--help")  Usage ;;
      *) Usage "Unknown option" $arg ;;
  esac
done

rc=${1:- tax_raw.rc}

if [ ! -f $rc ]
then
    echo "Please supply weights (.rc) file"
    exit 1
fi

scale=`kcalc $max_cost $min_cost -`

mx=`getcol -c 2 $rc | kstat -q -show max`
mn=`getcol -c 2 $rc | kstat -q -show min`

echo "# Created from $rc by normalise_tax_wgts.sh on `date`"
echo "# min_cost $min_cost max_cost $max_cost data min $mn data max $mx"
grep "^#" $rc
if [ $r_flag -eq 0 ]
then
    gawk -v min_cost=$min_cost -v scale=$scale -v mx="$mx" -v mn="$mn" '
    /^#/ {next;}
    $2 == "NA" {print $1, 2 * scale; next;}
    $2 < 0	{print; next}
            {print $1, min_cost + scale * ($2 - mn) / (mx - mn) }
' $rc
else
    # Reverse sense of entry
    gawk -v min_cost=$min_cost -v scale=$scale -v mx="$mx" -v mn="$mn" '
    $2 == "NA" {print $1, 2 * scale; next;}
    $2 < 0	{print; next}
            {print $1, min_cost + scale * (1.0 - ( ($2 - mn) / (mx - mn) ) ) }
' $rc
fi
