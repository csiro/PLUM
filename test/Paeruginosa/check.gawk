#!/bin/gawk
{   j = 2; 
    while (j < NF) {
        if ($j == "REACTION") {
            break;
        }
        j++;
    }
    middle = j;
    for (i = 1; i < middle; i++) {
        if ($i ~ /^[0-9.]*$/) {
            $i = $i * 1.0;
            $j = $j * 1.0;
        }
        if ($i != $j) {
            print;
            print i, j, $i, $j;
            break;
        }
        j++;
    }
}
