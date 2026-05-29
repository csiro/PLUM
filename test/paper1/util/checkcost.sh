#!/bin/bash
function run() {
    echo "$@"
    eval "$@"
}
Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo -n "Usage: `basename $0` [-lb #] [-ub #] [-sd fn] pld-fn" >& 2
    echo '
       See how how we can start cost
  Switches
       -lb: lower bound [1.0]
       -ub: upper bound [2000.0]
       -sd: sd file name [sd/glc_D.sd]
  Args
       pld-fn: the plum data file
' >&2
    exit 1
}

lb=1.0
ub=2000.0
sd_fn="sd/glc_D.sd"

while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-lb")  lb="$1" ; shift ;;
      "-ub")  ub="$1" ; shift ;;
      "-sd")  sd_fn="$1" ; shift ;;
      "-?")  Usage ;;
  "--help")  Usage ;;
      *) Usage "Unknown option" $arg ;;
  esac
done

if [ $# -lt 1 ]
then
    Usage "pld filename required"
fi
pld="$1"

dir="checkcost"
if [ ! -d $dir ]
then
    mkdir -v $dir
fi

res="checkcost.res"

plum_args="-v incr -rc tax.rc -sd $sd_fn"

cost=$lb

while true
do
    echo "Check cost $cost"
    diff=`kcalc $ub $lb - `
    termtitle "LB $lb UB $ub Cost $cost Diff $diff"
    
    sol="${dir}/${cost}.sol"
    rs="${dir}/${cost}.rs"
    model="${dir}/${cost}.pld"

    if [ ! -f $sol ]
    then
        run "plum $plum_args -o $sol -rs $rs -O $res -M $model -C $cost $pld"
    else
        echo "Already have sol for $cost"
    fi

    fp_cost=`kcalc 4 prec $cost`
    biomass=`grep "min_cost ${fp_cost}" $res | getcol biomass_flux`
    ibiomass=`kcalc $biomass rounddn`
    if [ "$biomass" == "" ]
    then
        echo "Can't find result for min_cost ${fp_cost}"
        break
    fi
    echo "Biomass flux is $biomass"

    if [ $ibiomass -lt 10 ]
    then
        echo "OK"
        lb=$cost
    else
        echo "Runaway"
        ub=$cost
    fi

    converge=0.01
    sgn=`kcalc $ub $lb - $converge - sgn`

    if [ $sgn -lt 0 ]
    then
        echo "Converged - LB $lb UB $ub"
        break
    fi
    cost=`kcalc $ub $lb - 2.0 / $lb +`
    echo "LB $lb UB $ub cost $cost"
done

echo "Done. Results in $res"
echo "min_cost biomass_flux"
getcol min_cost biomass_flux checkcost.res | sort -n
