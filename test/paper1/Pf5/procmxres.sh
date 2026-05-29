
res=${1:-"plummx-mr.res"}
echo "Using results file $res"

getcol react_cost_fn max_react_cost which_obj obj cost tau error $res | sort | sed -e 's|rc/tax_gfga-||' -e 's/\.rc//' | kstat -vn 3 -n 4 -H obj cost tau error

tmp="/tmp/procres.tmp"

getcol react_cost_fn max_react_cost which_obj cost seed $res | sort | sed -e 's|rc/tax_gfga-||' -e 's/\.rc//' | choosemin -n 3 | gawk '
/tax/ && $2 < 0 {print "out/out-tax.rc-" $4 ".out"; next}
/tax/ {print "out/out-tax.rc-" int($2) "-" $4 ".out"; next}
/-/   {print "out/out" $1 "-" $4 ".out"}' > $tmp

for f in `cat $tmp`
do
    echo -n "$f	"
    kjoin -c tax.rc $f | getcol 2 | sum
done
