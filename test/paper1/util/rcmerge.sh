#!/bin/bash

Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo -n "Usage: `basename $0` rc1 rc2" >& 2
    echo '
    Combine two reaction cost files
    Chooses the minimum cost from the two files
  Switches
      -r: reverse - choose maximum cost
     -v1: duplicate reactions in rc1, adding "__rev"
     -v2: duplicate reactions in rc2, adding "__rev"
  Args
      rc[12]: React cost files 
' >&2
    exit 1
}

choosemin_arg=""
minmax="min"
v1_flag=0
v2_flag=0
while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-r")  choosemin_arg="-r" ; minmax="max" ;;
      "-v1")  v1_flag=1 ;;
      "-v2")  v2_flag=1 ;;
      "-?")  Usage ;;
  "--help")  Usage ;;
      *) Usage "Unknown option" $arg ;;
  esac
done

if [ "$2" == "" ]
then
    Usage "2 react cost files required"
fi

rc1="$1"
rc2="$2"
tmp=/tmp/rcmerge.tmp

echo -n "" > $tmp
if [ $v1_flag -eq 1 ]
then
    gawk '!/^#/ {print $1 "__rev", $2}' $rc1 | grep -v '__rev__rev' >> $tmp
fi
if [ $v2_flag -eq 1 ]
then
    gawk '!/^#/ {print $1 "__rev", $2}' $rc2 | grep -v '__rev__rev' >> $tmp
fi
echo "# Merge (${minmax}) of $rc1 and $rc2 using rcmerge.sh on `date`"
echo "# ${rc1}:"
grep '^#' $rc1 
echo "# ${rc2}:"
grep '^#' $rc2
cat "$@" $tmp | sort | choosemin $choosemin_arg

