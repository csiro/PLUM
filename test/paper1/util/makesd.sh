#!/bin/bash

Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo -n "Usage: `basename $0` [-z #] [-t #] biolog.csv target-dir base.sd" >& 2
    echo '
    Create sd files from entries in biolog csv file
  Args
    biolog.csv: biolog CSV file
    target-dir: directory to put new sd files
       base.sd: Media sd appended to all sd files
  switches
        -z: no-growth value [0.2]
' >& 2
    exit 1
}

zero=0.2

while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-z")  zero="$1" ; shift ;;
      "-?")  Usage ;;
  "--help")  Usage ;;
      *) Usage "Unknown option" $arg ;;
  esac
done

if [ $# -lt 3 ]
then
    Usage "3 filenames required"
fi

ex_col=3
name_col=2
score_col=7
#ex_col=2
#name_col=1
#score_col=5

biolog="$1"
outdir="$2"
basesd="$3"

growth="growth.lis"
nogrowth="nogrowth.lis"

if [ ! -f "$biolog" ]
then
    echo "File \"${biolog}\" does not exist"
    exit 1
fi
if [ ! -d "$outdir" ]
then
    echo "Directory \"$outdir\" does not exist - creating"
    mkdir -v "$outdir"
fi

tmp="/tmp/makesd.tmp"

echo -n "" > $growth
echo -n "" > $nogrowth

echo "Using    ex col `head -1 "$biolog" | getcol -F, $ex_col`"
echo "Using  name col `head -1 "$biolog" | getcol -F, $name_col`"
echo "Using score col `head -1 "$biolog" | getcol -F, $score_col`"

gawk -F, -v ex_col=$ex_col -v name_col=$name_col '
NR > 1 && $ex_col != "NA"  && $1 != ""  {
       ex = $ex_col;
       name = $name_col;
       sub("^EX_","",ex);   
       sub("__","_",ex);   
       sub("_e$","__91__e__93__",ex);   
       sub("[(]","__91__",ex);
       sub("[)]","__93__",ex);
       gsub("\"","",name);
       print ex, "\"" name "\"" 
}' "$biolog" > $tmp

gawk -F, -v ex_col=$ex_col -v name_col=$name_col -v score_col=$score_col '
NR > 1 && $ex_col != "NA" && $1 != "" {
       ex = $ex_col;
       name = $name_col;
       score = $score_col;
       sub("^EX_","",ex);   
       sub("__","_",ex);   
       sub("_e$","__91__e__93__",ex);   
       id = gensub("[(].*","", 1, ex);
       sub("[(]","__91__",ex);
       sub("[)]","__93__",ex);
       gsub("\"","",name);
       print ex ";" id ";" score ";" name;
}' "$biolog" > $tmp

IFS=";"
while read -r line
do
    fld=( $line )
    ex=${fld[0]}
    id=${fld[1]}
    score=${fld[2]}
    name=${fld[3]}
    sdfn="${outdir}/${id}.sd"
    pos=`kcalc $score $zero - 1000 x int`
    echo "Name $name score is $score pos is $pos"
    if [ $pos -gt 0 ]
    then
        echo "  Growth"
        echo $sdfn >> $growth
    else
        echo "  Nogrowth"
        echo $sdfn >> $nogrowth
    fi
    #score=`kcalc $score $zero - 0 max`
    echo "Write $sdfn with $ex ( $score )"
    echo "biolog $score" > $sdfn
    echo "$ex 	0 10 true 	# ${name}" >> $sdfn
    cat $basesd >> $sdfn
done < $tmp
unset IFS

echo "Done. Created `wc -l < $tmp` sd files in directory ${outdir}"
echo "Growth sd's in $growth (`wc -l < $growth`); no-growth sd's in $nogrowth  (`wc -l < $nogrowth`)"

