#!/bin/bash
#!/bin/bash

Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo -n "Usage: `basename $0` [-v] file.pld" >& 2
    echo '
    Find EX and DM reactions in the PLD
    - EX reactions have no reactants
    - DM reactions have no products
  Switches
      -v: Reverse - select everything _except_ EX and DM    
  Args
      file.pld: plum data file
' >& 2
    exit 1
}

reverse=0
while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-v")  reverse=1 ;;
      "-?")  Usage ;;
  "--help")  Usage ;;
      *) Usage "Unknown option" $arg ;;
  esac
done

pld="$1"
if [ -z $pld ]
then
    Usage "PLD file required"
fi

gawk '
/^REACT/ {
         is_ex = 1;
         is_dm = 1;
         for (k = 6; k <= NF; k += 2) {
             if ($k < 0) {is_ex = 0;}
             if ($k > 0) {is_dm = 0;}
             if (k+1 <= NF) {
                k1 = k + 1
                if ($k1 == "__fullname__") {
                   break;
                }
             } 
         }
         if (is_ex == 1 || is_dm == 1) {
            if (!reverse) {
               print $0;  
            }
         }
         else {
            if (reverse) {
               print $0;
            }
         }
}' $pld
