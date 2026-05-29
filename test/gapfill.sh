#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --export=NONE

if [ "$1" == "" ]
then
    echo "Data filename required"
    exit 1;
fi

dat=$1
shift
base=`basename $dat .dat`

solarg="0"
tl="0"
abs="0"
fm="0"

for (( i=0; i<=$#; i++ ))
do
    arg="${!i}"
    case $arg in
        "-v")  $(( i++ )); solarg="${!i}" ;;
        "-t")  $(( i++ )); tl="${!i}" ;;
        "-a")  $(( i++ )); abs="${!i}" ;;
        "-fm")  $(( i++ )); fm="${!i}" ;;
    esac
done
solver="cts"
case $solarg in
    "0") solver="cts" ;;
    "1") solver="int" ;;
    "2") solver="cmb" ;;
esac
obj="rel"
case $abs in
    "0") obj="rel" ;;
    "1") obj="abs" ;;
    "f*") obj="rel" ;;
    "t*") obj="abs" ;;
esac
case $fm in
    "1") obj="${obj}-fm" ;;
    "t*") obj="${obj}-fm" ;;
esac

base="out/${base}-${solver}-${obj}-${tl}"

dig=${base}.dig
out=${base}.out

source ../module.sh

echo "Gapfill starting on `date`"
echo "Output base is $base"

cmd="../build/bin/gapfill -g $dig -o $out -O res -p 8 $dat $@"
echo $cmd
$cmd

echo "Gapfill finished on `date`"

