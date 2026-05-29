#!/bin/bash

# Display message and wait for a button-press
function pause() {
 read -s -n 1 -p "Press any key to continue . . ."
 echo ""
}
# Run a command
function run() {
    echo "$@"
    "$@"
    if [ $? -ne 0 ]
    then
        pause
    fi
}


Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo "Usage: `basename $0` [-c fn] datafn [growth-sd [nogrowth-sd]]" >& 2
    echo "    Runs multi-experiment solver" >& 2
    echo "  Args" >& 2
    echo "    datafn: plum format scenario" >& 2
    echo "    growth-sd: list of supply/demand files for growth media" >& 2
    echo "    nogrowth-sd: list of supply/demand files for zero-growth media" >& 2
    echo "  Switches" >& 2
    echo "     -c: Config filename. Put all plum switch values here" >& 2
    echo "  " >& 2
    echo "  - Run plum on scenario with multi-experiment INT solver" >& 2
    echo "    - Produce new model" >& 2
    echo "  - Run plum on new model for each sd in growth-sd and no-growth" >& 2
    echo "    - Check ranking" >& 2
    echo "  " >& 2
    exit 1
}

config="config"

while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-c")  config="$1" ; shift ;;
      "-?")  Usage ;;
  "--help")  Usage ;;
      *) Usage "Unknown option" $arg ;;
  esac
done

if [ $# -lt 1 ]
then
    Usage
fi

dat="$1"
growth_sd_lis="${2:-growth.lis}"
nogrowth_sd_lis="${3:-nogrowth.lis}"

# Other params
solver="int3"
out_dir="out"



if [ ! -f $config ]
then
    touch $config
fi

echo "          Dat is $dat"
echo "  Growth list is $growth_sd_lis"
echo "Nogrowth list is $nogrowth_sd_lis"


tmp="/tmp/runmulti1.tmp"
tmp2="/tmp/runmulti2.tmp"
tmp3="/tmp/runmulti3.tmp"
tmp4="/tmp/runmulti4.tmp"
tmp5="/tmp/runmulti5.tmp"
tmp6="/tmp/runmulti6.tmp"

base="`basename $dat .dat`"
multi_sol="${out_dir}/multi.sol"
multi_cts_sol="${out_dir}/multi-cts.sol"
newmodel="${out_dir}/multi.dat"

res_multi="${out_dir}/res.multi"
res_new="${out_dir}/res.new"
res_pa01="pa01/res"

if [ ! -d ${out_dir} ]
then
    mkdir -v ${out_dir}
fi

# Save config
cp -v $config ${out_dir}

if [ ! -f $newmodel ] || [ $dat -nt $newmodel ] || [ $config -nt $newmodel ]
then
    # Run plum on new scenario with multi-experiment INT solver
    colour YELLOW; echo "Run multi-INT solver on $dat"; colour RESET
    echo -n "" > $res_multi
    run plum -c $config data=$dat -v $solver -SD $growth_sd_lis -O $res_multi -o $multi_sol -oc $multi_cts_sol -M+ $newmodel
else
    echo "Multi solution is up to date"
fi

# Run plum on new model for each sd in growth-sd and nogrowth-sd 
colour YELLOW; echo "Run individual carbon sources for new model"; colour RESET
echo -n "" > $tmp2
for sd in `cat $growth_sd_lis $nogrowth_sd_lis`
do
    base="`basename $sd .sd`"
    sol="${out_dir}/${base}-new.sol"
    if [ $newmodel -nt $sol ]
    then 
        colour YELLOW; echo "  Run plum for $base (new model)"; colour RESET
        grep -v " sd_fn $sd " $res_new > $tmp
        mv $tmp $res_new
        run plum -q -c $config data=$newmodel -v cts -sd $sd -O $res_new -o $sol
    else
        echo "  Sol ($sol) for $base is up to date"
    fi
    
    grep -H biolog $sd | sed 's/:biolog/ /' >> $tmp2
done

# - Check ranking
colour YELLOW; echo "Check ranking"; colour RESET

dig="${out_dir}/multi.dig"
echo "H \"Correlation\"" 	> $dig
echo "T \"`date`\"" 		>> $dig
echo "X \"Biolog equiv\"" 	>> $dig
echo "Y \"Actual\"" 		>> $dig

# Find biomass flux of top biolog media
mx_sd=`sort -nr -k 2,2 $tmp2 | head -1 | getcol 1`
mx=`grep "sd_fn $mx_sd " $res_new | getcol biomass_flux`
biolog-equiv.sh $mx $growth_sd_lis $nogrowth_sd_lis > $tmp3
echo "L \"Correlation Line\"" 	>> $dig
echo "MD 0 0 $mx $mx" 		 	>> $dig
echo "" 		 				>> $dig

csv="${out_dir}/multi.csv"
echo "Time,`date`" > $csv
echo "Data,$dat" >> $csv
echo "Config" >> $csv
sed 's/\s\+/,/' config >> $csv

getcol sd_fn biomass_flux $res_new > $tmp
getcol sd_fn biomass_flux $res_pa01 > $tmp6

for lis in $growth_sd_lis # $nogrowth_sd_lis
do
    colour YELLOW; echo "$lis"; colour RESET
    
    kjoin $tmp3 $lis -o1 > $tmp5

    cat $tmp2 | kjoin - $lis -o1 | sort -k1,1 | sort -snr -k2,2 | kjoin - $tmp5 | kjoin - $tmp | kjoin - $tmp6 > $tmp4
    echo "C-Source Biolog BiologEquiv NewModel PA01"

    gawk '{printf ("%15s %4d %6.4f %6.4f %6.4f\n",
                      $1, $2,  $4,  $6,   $8); }' $tmp4
    
    echo "" >> $csv
    echo "$lis" >> $csv
    echo "Carbon Source,Biolog,BiologEquiv,NewModel,PA01" >> $csv
    gawk '{print $1 "," $2 "," $4 "," $6 "," $8 }' $tmp4 >> $csv

    echo ""
    echo "Kendall-tau for NewModel and PA01:"
    getcol 4 6 8 $tmp4 | kstat -n 3 -H Biolog NewModel PA01 -show tau
    
    echo "" >> $csv
    echo "Kendall-tau" >> $csv
    getcol 4 6 8 $tmp4 | kstat -n 3 -H Biolog NewModel PA01 -show tau | sed 's/ \{2,\}/,/g' >> $csv

    gawk '{print "M", $4, $6, "\nLP", $4, $6, $1}' $tmp4 >> $dig
done

echo "Done. Solution in $multi_sol. New model in $newmodel"
echo "Stats in $csv - view with"
echo "open $csv"
echo "Dig in $dig - view with"
echo "open $dig"
