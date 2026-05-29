#!/bin/bash

# Run a command
function run() {
    echo "$@"
    "$@"
}
# Display message and wait for a button-press
function pause() {
 read -s -n 1 -p "Press any key to continue . . ."
 echo ""
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
    echo "  - Run plum for each sd in nogrowth-sd" >& 2
    echo "  - Run plum for each sd in growth-sd" >& 2
    echo "  - Form union of all reactions used in growth-sd" >& 2
    echo "  - Create new scenario with selected reactions " >& 2
    echo "    - Reactions with cost <= 1.0 in datafn are preserved " >& 2
    echo "    - Other reactions come from union  " >& 2
    echo "  - Run plum on new scenario with multi-experiment INT solver" >& 2
    echo "    - Produce new model" >& 2
    echo "  - Run plum on new model for each sd in growth-sd and no-growth" >& 2
    echo "    - Check ranking" >& 2
    echo "  " >& 2
    exit 1
}

config="config"
if [ ! -f $config ]
then
    touch $config
fi

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


solver="CTS"

tmp="/tmp/runintmulti1.tmp"
tmp2="/tmp/runintmulti2.tmp"
tmp3="/tmp/runintmulti3.tmp"
tmp4="/tmp/runintmulti4.tmp"
tmp5="/tmp/runintmulti5.tmp"
tmp6="/tmp/runintmulti6.tmp"
tmp7="/tmp/runintmulti7.tmp"
tmp8="/tmp/runintmulti8.tmp"

base="`basename $dat .dat`"
union_dat="out/${base}-U.dat"
union_react="out/${base}-U.react"
multisol="out/multi.sol"
newmodel="out/multi.dat"

res_indiv="out/res.indiv"
res_multi="out/res.multi"
res_new="out/res.new"
res_union="out/res.union"
res_pa01="pa01/res"

if [ ! -d out ]
then
    mkdir -v out
fi

# Save config
cp -v $config out

# Run plum for each sd in no-growth-sd
colour YELLOW; echo "Run individual carbon sources"; colour RESET

touch $res_indiv
for sd in `cat $nogrowth_sd_lis`
do
    base="`basename $sd .dat`"
    sol="out/${base}-indiv.sol"
    # Run again if anthing changed
    if [ ! -f $sol ] || [ $dat -nt $sol ] || [ $sd -nt $sol ] || [ $config -nt $sol ] 
    then
        colour YELLOW; echo "  Run plum for $base"; colour RESET
        # Remove the result for this sd_fn from individual results
        grep -v " sd_fn $sd " $res_indiv > $tmp
        mv $tmp $res_indiv
        
        run plum -q -c $config data=$dat -sd $sd -O $res_indiv -o $sol -v $solver
    else 
        colour YELLOW; echo "  Have indiv sol for $base"; colour RESET; 
    fi
done

sol_list=""
changed=0
# Run plum for each sd in growth-sd
for sd in `cat $growth_sd_lis`
do
    base="`basename $sd .dat`"
    sol="out/${base}-indiv.sol"
    sol_list="$sol_list $sol"
    # Run again if anthing changed
    if [ ! -f $sol ] || [ $dat -nt $sol ] || [ $sd -nt $sol ] || [ $config -nt $sol ] 
    then
        changed=1
        colour YELLOW; echo "  Run plum for $base"; colour RESET
        # Remove old result for sd
        grep -v " sd_fn $sd " $res_indiv > $tmp
        mv $tmp $res_indiv
        run plum -q -c $config data=$dat -sd $sd -O $res_indiv -o $sol -v $solver
    else 
        colour YELLOW; echo "  Have indiv sol for $base"; colour RESET; 
    fi
done

# Save biomass fluxes for each sd
getcol sd_fn biomass_flux $res_indiv > $tmp2

# Form union of all reactions used by growth 
if [ $changed -eq 1 ] || [ ! -f $union_dat ]
then
    colour YELLOW; echo "Form union of reactions"; colour RESET
    echo "sol list is $sol_list" 
    cat $sol_list | grep -v '^#' | getcol 1 | sort | uniq > $union_react

    # Create new scenario with selected reactions
    echo "# Created by `basename $0` on `date` from" > $union_dat
    echo "#   datafile: $dat" 						>> $union_dat
    echo "#   sd files: `unwrap $growth_sd_lis`" 	>> $union_dat
    grep ^MET $dat									>> $union_dat
    # - Reactions with cost <= 1.0 in datafn are preserved 
    gawk '/^\s*REACT/ && $3 <= 1' $dat 				>> $union_dat
    # - Other reactions come from union
    # Copy reactions from dat if they are in union, and cost more than 1
    gawk '/^REACT/ && $3 > 1 {$3 = 10; print}' $dat | kjoin $union_react - -j2 2 -o2 >> $union_dat

    #ause
fi

# Run plum on new scenario with multi-experiment INT solver
colour YELLOW; echo "Run multi-INT solver on $union_dat"; colour RESET
if [ ! -f $res_multi ] || [ ! -f $multisol ] || [ ! -f $newmodel ] || [ $union_dat -nt $newmodel ]
then
    echo -n "" > $res_multi
    run plum -c $config data=$union_dat -v INT -x -1 -SD $growth_sd_lis -bf $tmp2 -O $res_multi -o $multisol -M+ $newmodel
else 
    colour YELLOW; echo "  Have sol for $union_dat"; colour RESET; 
fi

# Run plum on new model for each sd in growth-sd and nogrowth-sd 
colour YELLOW; echo "Run individual carbon sources for new model"; colour RESET
echo -n "" > $tmp3
for sd in `cat $growth_sd_lis $nogrowth_sd_lis`
do
    base="`basename $sd .dat`"
    sol="out/${base}-new.sol"
    md="out/${base}-new.md"
    mds="out/${base}-new.mds"
    if [ ! -f $sol ] || [ $newmodel -nt $sol ] || [ $sd -nt $sol ] || [ $config -nt $sol ] 
    then
        colour YELLOW; echo "  Run plum for $base (new model)"; colour RESET
        grep -v " sd_fn $sd " $res_new > $tmp
        mv $tmp $res_new
        run plum -q -c $config data=$newmodel -sd $sd -O $res_new -o $sol -md $md
        # Created sorted biomass duals file
        gawk '$2 == 1' $md | sort -gr -k3,3 > $mds
    else 
        colour YELLOW; echo "  Have new sol for $base"; colour RESET; 
    fi

    # Now run plum on union data to check flux.
    sol="out/${base}-union.sol"
    if [ ! -f $sol ] || [ $union_dat -nt $sol ] || [ $sd -nt $sol ] || [ $config -nt $sol ] 
    then
        colour YELLOW; echo "  Run plum for $base (union data)"; colour RESET
        grep -v " sd_fn $sd " $res_union > $tmp
        mv $tmp $res_union
        run plum -q -c $config data=$union_dat -sd $sd -O $res_union -o $sol
    else 
        colour YELLOW; echo "  Have union sol for $base"; colour RESET; 
    fi

    grep -H biolog $sd | sed 's/:biolog/ /' >> $tmp3
done

# - Check ranking
colour YELLOW; echo "Check ranking"; colour RESET

# Find biomass flux of top biolog media
mx_sd=`sort -nr -k 2,2 $tmp3 | head -1 | getcol 1`
mx=`grep "sd_fn $mx_sd " $res_new | getcol biomass_flux`
biolog-equiv.sh $mx $growth_sd_lis $nogrowth_sd_lis > $tmp5

csv="out/multi.csv"
echo "Time,`date`" > $csv
echo "Data,$dat" >> $csv
echo "Config" >> $csv
sed 's/\s\+/,/' config >> $csv

getcol sd_fn biomass_flux $res_new > $tmp
getcol sd_fn biomass_flux $res_union > $tmp4
getcol sd_fn biomass_flux $res_pa01 > $tmp8

for lis in $growth_sd_lis $nogrowth_sd_lis
do
    colour YELLOW; echo "$lis"; colour RESET
    
    kjoin $tmp5 $lis -o1 > $tmp7

    cat $tmp3 | kjoin - $lis -o1 | sort -k1,1 | sort -snr -k2,2 | kjoin - $tmp7 | kjoin - $tmp | kjoin - $tmp2 | kjoin - $tmp4 | kjoin - $tmp8 > $tmp6
    echo "C-Source Biolog BiologEquiv NewModel Indiv Union PA01"

    gawk '{printf ("%15s %4d %6.4f %6.4f %6.4f %6.4f %6.4f\n",
                      $1, $2,  $4,  $6,   $8,   $10,  $12); }' $tmp6
    
    echo "" >> $csv
    echo "$lis" >> $csv
    echo "Carbon Source,Biolog,BiologEquiv,NewModel,Indiv,Union,PA01" >> $csv
    gawk '{print $1 "," $2 "," $4 "," $6 "," $8 "," $10 "," $12}' $tmp6 >> $csv

    echo ""
    echo "Kendall-tau for NewModel and PA01:"
    getcol 4 6 12 $tmp6 | kstat -n 3 -H Biolog NewModel PA01 -show tau
    
    echo "" >> $csv
    echo "Kendall-tau" >> $csv
    getcol 4 6 12 $tmp6 | kstat -n 3 -H Biolog NewModel PA01 -show tau | sed 's/ \{2,\}/,/g' >> $csv
done

echo "Done. Solution in $multisol. New model in $newmodel"
echo "Stats in $csv - view with"
echo "open $csv"
