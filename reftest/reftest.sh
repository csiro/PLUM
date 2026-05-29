#!/bin/bash

Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo "Usage: `basename $0` [-x] [reftest.sh]" >& 2
    echo "       Run tests " >&2
    echo "       Input file contains bash commands to do runs "
    echo "          and diff commands to test output" >&2
    echo "  Switches" >&2
    echo "    -x: Output an example reftest.sh" >&2
    echo "  Args" >&2
    echo "    testfile: RUN and TEST commands [reftest.dat]" >&2
    exit 1
}

Example () {
    cat << __EOF__
# Prep
if [ ! -d out ]; then mkdir out; else rm -v out/*; fi

# Test 1 - CTS solver
plum -s 1 -O out/res data/PA01.dat -sd data/sd-4abut.dat -o out/test1.out
diff out/test1.out ref/test1.out

# Test 2 - INT solver
plum -s 1 -O out/res data/PA01.dat -sd data/sd-4abut.dat -o out/test2.out -v INT
diff out/test2.out ref/test2.out
__EOF__
    exit 1
}

function run() {
    echo "$@"
    "$@"
}

# Process switches
while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-x")  Example ;;
      "-?")  Usage ;;
  "--help")  Usage ;;
      *) Usage "Unknown option" $arg ;;
  esac
done

testfile="${1:-reftest.dat}"
if [ ! -r $testfile ]
then
    echo "Couldn't open file $testfile" >&2
    exit
fi
    
out="reftest.out"
echo "reftest run on `date`" > $out

while read -r line
do
    if [ "$line" == "" ]
    then
        echo "================================================" >> $out
    else
        colour YELLOW; echo $line; colour RESET
        echo "------------------------------------------------" >> $out
        echo $line >> $out
        eval "$line" | tee -a $out
    fi
done < $testfile

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 
colour YELLOW
echo "Done. Output in $out. View with"
colour RESET
echo "less $out"
