#!/bin/bash

res=${1:-"res"}

echo "Key react_cost_fn which_obj "
getcol -c react_cost_fn which_obj obj cost tau error $res | sort | kstat -vn 2 -n 4 -H obj cost tau error

declare -a key
tmp=/tmp/procres.tmp
getcol -c react_cost_fn which_obj $res | sort | uniq | while IFS='' read -r line
do
    key=( $line )
    #echo "line is $line key is (${key[0]}, ${key[1]})"
    echo "Displayed Pareto is cost tau error react_cost_fn which_obj seed"
    grep " react_cost_fn ${key[0]}" $res | grep " which_obj ${key[1]} " | getcol -c cost tau error react_cost_fn which_obj seed - | sort | pareto -m2 
done

echo "Overall cost/tau pareto front is "
echo "Displayed Pareto is cost tau error react_cost_fn which_obj seed"
getcol cost tau error react_cost_fn which_obj seed $res | sort -n | pareto -m2 | cat -n
echo ""
echo "Overall cost/error pareto front is "
echo "Displayed Pareto is cost error tau react_cost_fn which_obj seed"
getcol cost error tau react_cost_fn which_obj seed $res | sort -n | pareto | cat -n
