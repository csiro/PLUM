#!/bin/bash
function run() {
    echo "$@"
    eval "$@"
}
# Display message and ask for a reply with return. Return keypress in $ask
ask=""
function ask_wait() {
 read -p "$1 " ask
 echo ""
}
Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo -n "Usage: `basename $0` [-a str] [-d dtr] pld-fn sd-lis" >& 2
    echo '
       Run all sd files in sd-lis
       Config defined by "testall.cfg"
  Switches
       -a: extra arg passed to plum. Can be repeated
       -d: output directory [testall]
       -O: res file (appended) [none]
  Args
       pld-fn: the plum data file
       sd-lis: list of sd-files to use
' >&2
    exit 1
}

plum_args="-c testall.cfg"
dir="testall"
myres="none"

while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-a")  plum_args="$plum_args $1" ; shift ;;
      "-d")  dir="$1" ; shift ;;
      "-O")  myres="$1" ; shift ;;
      "-?")  Usage ;;
  "--help")  Usage ;;
      *) Usage "Unknown option" $arg ;;
  esac
done

if [ $# -lt 2 ]
then
    Usage "2 filenames required"
fi

pld="$1"
sd_lis="$2"

res="${dir}/testall.res"
union="${dir}/union.txt"

if [ ! -d $dir ]
then
    mkdir -v $dir
fi

if [ -f $res ]
then
    ask_wait "Remove old results? "
    if [ "$ask" == "y" ] || [ "$ask" == "Y" ]
    then
        rm -v $res
        rm -rv $dir
        mkdir -v $dir
    fi
fi

biolog="${dir}/biolog.lis"
echo -n "" > $biolog

all_sols=""

for sd in `cat $sd_lis`
do
    if [[ "$sd" == "" ]] || [[ "$sd" =~ ^"#" ]]
    then
        continue
    fi
    b=`basename $sd .sd`
    sol="${dir}/${b}.sol"
    termtitle "sd $sd"

    echo -n "$sd " >> $biolog
    grep biolog $sd >> $biolog

    if [ ! -f $sol ]
    then
        run plum $plum_args -sd $sd -o $sol -O $res $pld
    else
        echo "Already have sol for $sd"
    fi
    all_sols="$all_sols $sol"
done

cat $all_sols | getcol -c 1 3 | sort | uniq -c | getcol 2 1 3 - > $union

summary="${dir}/summary.txt"
echo "sd_fn biomass_flux biolog"
echo "sd_fn biomass_flux biolog" > $summary
getcol sd_fn biomass_flux $res | kjoin - $biolog | getcol 1 2 5 | sort -k3,3nr | tee -a $summary

tmp="/tmp/testall.tmp"
gawk '
BEGIN 	{num_match = 0; num_true_pos = 0; num_mismatch = 0; num_false_neg = 0; runaway = 0;}
NR > 1  { if ($3 <= 0.3) { 
             if ($2 > 0) {
                num_mismatch++;
             }
             else {
                num_match++;
             }
          }
          else {
             if ($2 > 0) {
                if ($2 > 10) {
                   runaway++;
                }
                else {
                   num_true_pos++;
                   num_match++;
                }   
             }
             else {
                num_mismatch++;
                num_false_neg++;
             }
          }
        }
END     {
            print ( \
                   "Growth match", num_match, \
                   "including", num_true_pos, "true positives", \
                   "plus", runaway, "run-away" \
            );
            print ( \
                   "Growth mismatch", num_mismatch,\
                   "including", num_false_neg, "false negatives" \
            );
        }
' $summary | tee $tmp

tau=`grep ^sd $summary | getcol 3 2 | kstat -n 2 -show tau -q | getcol 2`

if [ $myres != "none" ]
then
    match=`grep "Growth match" $tmp | getcol match -`
    true_pos=`grep "Growth match" $tmp | getcol including -`
    mismatch=`grep "Growth mismatch" $tmp | getcol mismatch -`
    false_neg=`grep "Growth mismatch" $tmp | getcol including -`

    echo "pld $pld sd_lis $sd_lis args $plum_args match $match mismatch $mismatch true_pos $true_pos false_neg $false_neg tau $tau" >> $myres
fi

cat $tmp >> $summary
echo "Kendall tau is $tau" | tee -a $summary

echo "Done. Results in $res. Union solution in $union. Summary in $summary"

