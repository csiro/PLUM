#!/bin/gawk
{   mult = -1; 
    printf ("REACTION R_%s 1", $3)
    for (i = 4; i <= NF; i++) {
        if ($i == "<=>" || $i == "-->") {
            mult = 1.0;
        }
        else if ($i == "+") {
            # Do nothing
        }
        else {
            fact = 1.0;
            if ($i ~ /^[0-9.]*$/) {
                fact = $i;
                i++;
            }
            fact = mult * fact;
            printf (" M_%s %f", $i, fact);
        }
    }
    printf ("\n");
}
            
           
