#!/bin/gawk
$1 ~ /^EX_/	{   
    gsub (":", "", $1);
    gsub ("EX_", "", $1);
    print "M_" $1, 0, $2 * -1.0;
}
