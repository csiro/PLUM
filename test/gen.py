#!/usr/bin/env python3
"""
Test data generator for metabolic gap-filling solver.

This module generates synthetic metabolic network data for testing flux balance
analysis (FBA) gap-filling algorithms. It creates random metabolites and reactions
with configurable parameters to simulate metabolic networks of varying complexity.

Examples
--------
Generate a network with 50 metabolites and 100 reactions:
    $ python gen.py 50 100

Generate with custom max metabolites per reaction and objective coefficient:
    $ python gen.py -m 8 -o 10 50 100

Notes
-----
- Metabolites are categorized into three groups: supply (20%), residual (20%),
  and intermediate (60%)
- Reactions are categorized: preferred reactions (20%) and regular reactions (80%)
- Stoichiometric coefficients are generated in quarters (0.25 increments)
- Each reaction is guaranteed to have both positive and negative coefficients
"""
import sys
import os
import fileinput
import random
from datetime import datetime
from itertools import chain

def die (message):
    """
    Print error message to stderr and exit the program.

    Parameters
    ----------
    message : str
        Error message to display before terminating.

    Raises
    ------
    SystemExit
        Always exits with code 1.

    Examples
    --------
    >>> die("Invalid input file")
    Invalid input file
    """
    print(message, file=sys.stderr)
    exit (1)

def usage (str = ""):
    """
    Display usage information and exit the program.

    Parameters
    ----------
    str : str, optional
        Additional message to display before usage information (default is "").

    Raises
    ------
    SystemExit
        Always exits with code 1.

    Notes
    -----
    Prints command-line syntax and argument descriptions for the test data
    generator, including:
    - n: number of metabolites in the network
    - m: number of reactions in the network
    - -m: maximum metabolites per reaction (default: 5)
    - -o: objective coefficient for preferred reactions (default: 0)
    """
    print (str)
    print (sys.argv[0], end=' ')
    print (''' [-m #] n m
        Generate test data for gap filling solver
    Args
        n: number of metabolites
        m: number of reactions
    Switches
        -m: max metabolites in a reaction (min 2) [5]
        -o: obj coeff for preferred reactions [0]
''')
    exit (1)


def main():
    """
    Main function to generate synthetic metabolic network test data.

    Parses command-line arguments and generates a random metabolic network with
    specified numbers of metabolites and reactions. The output is formatted for
    use with gap-filling solvers.

    Returns
    -------
    None
        Prints formatted network data to stdout.

    Raises
    ------
    SystemExit
        Exits with code 1 if invalid arguments are provided.

    Notes
    -----
    Metabolite categories:
    - First 20%: Supply metabolites with random supply values (2-10)
    - Next 20%: Residual metabolites with random residual values (2-10)
    - Remaining 60%: Intermediate metabolites

    Reaction categories:
    - First 20%: Preferred reactions with user-specified objective coefficient
    - Remaining 80%: Regular reactions with random objective coefficients (100-300)

    Stoichiometric coefficients:
    - Generated in increments of 0.25 (quarters)
    - Range: -2.0 to 1.75 (excluding 0)
    - Each reaction guaranteed to have mixed signs (both reactants and products)

    Examples
    --------
    Generate network with 50 metabolites and 100 reactions:
    >>> main()  # with sys.argv = ['gen.py', '50', '100']

    Generate with custom parameters:
    >>> main()  # with sys.argv = ['gen.py', '-m', '8', '-o', '10', '50', '100']
    """
    # parse argv
    n = 0
    m = 0
    max_coeff = 5
    pref_obj = 0
    while len(sys.argv) > 1:
        arg = sys.argv.pop(1)
        if arg[0] == '-':
            if arg == '-m' and len(sys.argv) > 1:
                max_coeff = int(sys.argv.pop(1))
            elif arg == '-o' and len(sys.argv) > 1:
                pref_obj = int(sys.argv.pop(1))
            else:
                usage("Unrecognised arg: " + arg)
        else:
            # Alternatively
            if n == 0:
                n = int(arg)
            elif m == 0:
                m = int(arg)
            else:
                usage("Too many args: " + arg)

    if (n == 0 or m == 0):
        usage()
    
    random.seed()
    now = datetime.now()
    
    print ("# Created by gen.py on", now.strftime("%d/%m/%Y %H:%M:%S"))
    print ("# num metabolites = " + str(n) + 
           ", num reactions = " + str(m) + 
           ", max coeff = " + str(max_coeff))
    max_coeff = 5
    pref_obj = 0
    # Write metabolites
    for k in range(n):
        # 1/5 are supply, 1/5 are residual, and rest are used in reactions
        supply = 0
        residual = 0
        if k < n / 5:
            # It is supplied 
            supply = random.randint (2, 10)
        elif k < 2 * n / 5:
            # It has resual
            residual = random.randint (2, 10)
            
        print ("MET met" + str(k+1), supply, residual)
        
    print ()
    # Write reactions with a random number of metabolites
    for k in range(m):
        # First 1/5 are 'known' reactions with 'preferred' obj value
        obj_coeff = pref_obj
        if k > m/5:
            obj_coeff = 100 * random.randint (1, 3)
        
        print ("REACTION react" + str(k+1), obj_coeff, end=" ")
        
        num_coeff = random.randint (2,max_coeff)
        all_pos = True
        all_neg = True
        
        choices = [*range(n)]
        
        for i in range(num_coeff):
            # Coeffs are in quarters, with random sign
            met = random.choice (choices) 
            choices.remove(met)
            met_coeff = random.randint (-8, 7)
            if (met_coeff == 0): # Exclude 0
                met_coeff = met_coeff + 1
            met_coeff = met_coeff / 4.0
            if met_coeff > 0:
                all_pos = False
            if met_coeff < 0:
                all_neg = False
            # On last one, if all are the same sign, swap it
            if i == (num_coeff - 1) and (all_neg or all_pos):
                met_coeff = -met_coeff
            
            print (" met"+str(met+1), met_coeff, end="")
        print ()
            
        

if __name__ == '__main__':
    main()
