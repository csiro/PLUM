#!/usr/bin/env python3
import sys
import os
import fileinput
from datetime import datetime # For date formatting
import csv

obj_value = 2000
default_out_fn = "out.pld"

include_ex = False
include_dm = False

def die (message):
    print (message, file=sys.stderr)
    exit (1)

def usage (msg = ""):
    global obj_value
    
    print (msg)
    print (sys.argv[0], end=' ')
    print ('''[-v #] react.csv met.csv [out.pld]
        Convert csv files into a PLD file. 
        (removes biomass reactions)
    Args
        react.csv: CSV file with reactions 
          met.csv: CSV file with metabolites
          out.pld: output file [''' + default_out_fn + ''']
    Switches
        -v: obj value of reactions [''' + str(obj_value) + ''']
       +ex: include exchange reactions (no reactants) [don't]
       +dm: include demand reactions (no products) [don't]
''')
    exit (1)
    
def fix_name (name):
    name = name.replace ("[", "__91__")
    name = name.replace ("]", "__93__")
    name = name.replace ("(", "__40__")
    name = name.replace (")", "__41__")
    name = name.replace ("-", "__45__")
    return name

def is_float (string):
    try:
        v = float(string)
        return True
    except ValueError:
        return False
    


def main():
    global obj_value
    global include_ex
    global include_dm
    
    # parse argv
    react_csv = ""
    met_csv = ""
    pld_fn = ""
    
    react_name_col = 0
    react_id_col = 1
    stoich_col = 2
    react_cols = max(react_name_col,react_id_col,stoich_col)+1

    met_name_col = 1
    met_id_col = 2
    met_formula_col = 3
    met_cols = max(met_name_col,met_id_col,met_formula_col)+1

    while len(sys.argv) > 1:
        arg = sys.argv.pop(1)
        if arg[0] == '-' or arg[0] == '+':
            if arg == '-v' and len(sys.argv) > 1:
                obj_value = float(sys.argv.pop(1))
            elif arg == '+ex':
                include_ex = True
            elif arg == '+dm':
                include_dm = True
            else:
                usage("Unrecognised arg: " + arg)
        else:
            # Alternatively
            if react_csv == "":
                react_csv = arg
            elif met_csv == "":
                met_csv = arg
            elif pld_fn == "":
                pld_fn = arg
            else:
                usage("Too many file names: " + arg)

    if met_csv == "":
        usage("Need filenames")
        
    try:
        react = open (react_csv, 'r')
    except IOError:
        die ('Can\'t open file ' + react_csv + ' for reading')
    try:
        met = open (met_csv, 'r')
    except IOError:
        die ('Can\'t open file ' + met_csv + ' for reading')

    if pld_fn == "":
        pld_fn = default_out_fn
    try:
        pld = open (pld_fn, 'w')
    except IOError:
        die ('Can\'t open file ' + pld_fn + ' for writing')
        
    dt_string = datetime.now().strftime("%-d %b %Y %H:%M:%S") 
    print ("# Created by csv2pld.py on " + dt_string, file=pld)
    print ("# from " + react_csv + ", " + met_csv, file=pld)
    print ("# incl. ex " + str(include_ex) + 
           " incl. dm " + str(include_dm), file=pld)

    react_reader = csv.reader (react)
    fld = next(react_reader)

    print ("Reaction name from " + fld[react_name_col])
    print ("Reaction   id from " + fld[react_id_col])
    print ("Stoichiometry from " + fld[stoich_col])
    print ("")
    
    met_reader = csv.reader (met)
    fld = next(met_reader)
    
    print ("Met    name from " + fld[met_name_col])
    print ("Met      id from " + fld[met_id_col])
    print ("Met formula from " + fld[met_formula_col])
    print ("")
    
    react_lines = []
    used_mets = dict()
    skipped_biomass = 0
    skipped_ex = 0
    skipped_dm = 0

    for row in react_reader:
        if len(row) < react_cols:
            print ("Ignore " + str(row))
            continue
        name = row[react_name_col]
        ident = fix_name (row[react_id_col])
        stoich = row[stoich_col]
        
        fld = stoich.split ()
        pairs = []
        reactant = True
        direction = 0
        num_reactants = 0
        num_products = 0

        mult = 1.0
        for k in range(len(fld)):
            if is_float(fld[k]):
                mult = float(fld[k])
            elif fld[k] == '-->':
                direction = 1
                reactant = False
            elif fld[k] == '<--':
                direction = -1
                reactant = False
            elif fld[k] == '<=>':
                direction = 2
                reactant = False
            elif fld[k] == '+':
                pass 
            else:
                met = fix_name (fld[k])
                if (reactant):
                    mult *= -1
                    num_reactants += 1
                else:
                    num_products += 1
                pairs.append ([mult, met])
                mult = 1.0
                used_mets[met] = 1

        if direction == 0:
            print ("No direction for reaction " + name)
            exit (2)
            
        if ident.startswith("biomass"):
            skipped_biomass += 1
        else:
            if direction == 1 or direction == 2:
                if not include_ex and num_reactants == 0:
                    skipped_ex += 1
                elif not include_dm and num_products == 0:
                    skipped_dm += 1
                else:
                    formula = ""
                    for mult,met in pairs:
                        formula += " " + met + " " + str(mult)
                        
                    react_lines.append ( \
                        "REACTION " + ident + " " + str(obj_value) + 
                        " 1000" + 
                        formula + " __fullname__ \"" + 
                        name + "\"")
                    
            if direction == -1 or direction == 2:
                # Include reverse flow
                if not include_ex and num_products == 0:
                    skipped_ex += 1
                elif not include_dm and num_reactants == 0:
                    skipped_dm += 1
                else:
                    formula = ""
                    for mult,met in pairs:
                        formula += " " + met + " " + str(-mult)
            
                    ident += "__rev"
                    react_lines.append ( \
                        "REACTION " + ident + " " + str(obj_value) + 
                        " 1000" + formula + " __fullname__ \"" + 
                        name + " (rev)\"")

    mets = dict()

    met_lines = []
    for row in met_reader:
        if len(row) < met_cols:
            print ("Ignore " + str(row))
            continue
        
        name = row[met_name_col]
        ident = fix_name (row[met_id_col])
        formula = row[met_formula_col]
        
        if (ident in used_mets):
            met_lines.append ("MET " + ident + " \"" + name + 
                              " (" + formula + ")\"")
        
    for line in met_lines:
        print (line, file=pld)
    for line in react_lines:
        print (line, file=pld)
        
    print ("Done. Output in " + pld_fn);
    print ("  " + str(len(met_lines)) + " metabolites");
    print ("  " + str(len(react_lines)) + " reactions");
    print ("  " + str(skipped_biomass) + " biomass reactions skipped");
    print ("  " + str(skipped_ex) + " EX reactions skipped");
    print ("  " + str(skipped_dm) + " DM reactions skipped");
        
if __name__ == '__main__':
    main()
    
