#!/bin/bash

Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo -n "Usage: `basename $0` [-x #] [-rc1 fn] pdl1-fn pdl2-fn" >& 2
    echo '
  Check how many of the reactions are common between pdl1 and pdl2

  Switches
        -x: threshold for "expensive" reaction [200]    
      -rc1: filename for react costs for pld1 [none]
      -rc2: filename for react costs for pld2 [none]
  Args
       pdl[12]: Name of plum-format scenario file"
' >&2
    exit 1
}

thresh=200
rc1="none"
rc2="none"
# Another opt method
while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-x")  thresh="$1" ; shift ;;
      "-rc1")  rc1="$1" ; shift ;;
      "-rc2")  rc2="$1" ; shift ;;
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

pld1="$1"
pld2="$2"

if [ ! -r $pld1 ]
then
    Usage "Couldn't open file $pld1"
fi
if [ ! -r $pld2 ]
then
    Usage "Couldn't open file $pld2"
fi

tmp1=/tmp/checkoverlap1.tmp
tmp2=/tmp/checkoverlap2.tmp
tmp3=/tmp/checkoverlap3.tmp
tmp4=/tmp/checkoverlap4.tmp
tmp5=/tmp/checkoverlap5.tmp

grep ^REACTION $pld1 | gawk '$3 > 1 {print $2, $3}'  > $tmp1
grep ^REACTION $pld2 | gawk '$3 > 1 {print $2, $3}'  > $tmp2

if [ "$rc1" != "none" ]
then
    kjoin $rc1 $tmp1 -a2 | getcol 1 2 > $tmp3
    mv $tmp3 $tmp1
fi
if [ "$rc2" != "none" ]
then
    kjoin $rc2 $tmp2 -a2 | getcol 1 2 > $tmp3
    mv $tmp3 $tmp2
fi

kjoin $tmp1 $tmp2 > $tmp3

count1=`wc -l < $tmp1`
count2=`wc -l < $tmp2`
count3=`wc -l < $tmp3`

gawk "\$2 > $thresh" $tmp1 | sort -k2,2nr -k1,1 > $tmp4
gawk "\$2 > $thresh" $tmp2 | sort -k2,2nr -k1,1 > $tmp5

count4=`wc -l < $tmp4`
count5=`wc -l < $tmp5`

echo "$count1 non-gene-assoc reactions in $pld1"
echo "$count2 non-gene-assoc reactions in $pld2"
echo "$count3 reactions in common"
echo ""
echo "$count4 reactions in $pld1 with cost > $thresh"
echo "$count5 reactions in $pld2 with cost > $thresh"
echo ""
echo "Expensive reactions ($tmp4, $tmp5):"
paste $tmp4 $tmp5


