
dir=${1:-"plummx"}
ref=${2:-"PA01.pld"}
ref_tau=${3:- "-0.0435"}

res="$dir/res"
echo "Using results file $res"

filelist="`getcol data $res | sort | uniq | myfold`"
filelist="$filelist `getcol react_cost_fn $res | sort | uniq | myfold`"
filelist="$filelist $res"

tmp="/tmp/procmxres.tmp"

echo "Pareto front"
echo "Over all"
echo "  Cost Tau file"
grep '^# seed ' $dir/descr-*.out | sed 's/:/ /' | getcol cost tau 1 | pareto -m2

for rc in `getcol react_cost_fn $res | sort | uniq`
do
    rc=`basename $rc .rc`
    echo ""
    echo "Over $rc"
    echo "  Cost Tau file"
    grep '^# seed ' $dir/descr-*.out | grep "$rc\.rc" | sed 's/:/ /' | getcol cost tau 1 | pareto -m2 | tee $tmp
    echo ""
    for f in `getcol 3 $tmp`
    do
        echo 
        echo "$f "
        filelist="$filelist $f"
        echo "  Cost profile (to nearest 100)"
        echo -n "   "
        grep "full name" < $f | getcol cost |  kcalc =1 100 roundn | sort -n |  uniq -c | myfold
        echo "  Top 10 cost reactions"
        grep "full name" < $f | sort -k3,3n | tail -10 | gawk '{print "  " $0}'
    done
done

dig="$dir/procmxres.dig"
echo "T \"Produced by procmxres.sh from dir $dir on `date`\"" > $dig
echo "X Cost" >> $dig
echo "Y Tau" >> $dig
for rc in `getcol react_cost_fn $res | sort | uniq`
do
    rcb=`basename $rc .rc`
    echo "H \"$rcb\"" >> $dig
    
    if [ "$ref" != "" ] && [ -f $ref ]
    then
        ref_cost="`grep "^REACT" $ref | getcol 2 | kjoin - $rc | getcol 3 | sum`"
        echo "Ref cost from $ref using $rcb is $ref_cost"
        echo "" >> $dig
        echo "L \"Ref cost from $ref\"" >> $dig
        echo "M $ref_cost $ref_tau" >> $dig
    fi

    for obj in `getcol which_obj $res | sort -n | uniq`
    do
        name=`objname.sh $obj`
        echo "" >> $dig
        echo "L \"Obj $obj: $name\"" >> $dig

        grep '^# seed ' $dir/descr-*.out | grep $rc | grep "which_obj $obj " | sed 's/:/ /' | getcol cost tau 1 | pareto -m2 > $tmp
        gawk '{print "LP", $1, $2, "\"", $3  "\"\n" "D", $1, $2}' $tmp >> $dig

        interesting="`getcol 3 $tmp | myfold`"
        filelist="$filelist $interesting"
    done
    echo "W" >> $dig
    echo "WIPE" >> $dig
done

zip="$dir/plummx.zip"
echo ""
echo "filelist is "
echo "  $filelist"
zip -q $zip $filelist
echo "Interesting files in $zip"

echo "Also check out dig file $dig with"
echo "dig $dig &"
