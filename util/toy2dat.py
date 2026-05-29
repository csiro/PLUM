#!/usr/bin/env python3

import sys
import os
import fileinput
import openpyxl

std_names = ["EP", "PC", "CE"]

def die (message):
    print (message, file=sys.stderr)
    exit(1)

def usage (str = ""):
    print (str)
    print (sys.argv[0], end=' ')
    print ('''[-f #] in.xlsx 
        Convert a toy prob xlsx into gapfill data
        Writes to stdout
    Args
        file to process
    Switches
        -f: column for flux (1-based) [6]
''')
    exit(1)

def main ():
    global std_names
    flux_col = 5
    dir_col = 4
    filename = ""
    while len(sys.argv) > 1:
        arg = sys.argv.pop(1)
        if (arg[0] == '-'):
            if arg == '-f' and len(sys.argv) > 1:
                flux_col = int(sys.argv.pop(1)) - 1
            else:
                usage("Unrecognised arg: " + arg)
        else:
            # Save for file processing
            if filename == "":
                filename = arg
            else:
                usage ("Too many filenames")
                
    if filename == "":
        usage ("Filename required")

    to_compart_col = flux_col + 1
    from_compart_col = to_compart_col + 1
    dir_col = flux_col - 1
    
    wb = openpyxl.load_workbook(filename = filename)
    ws = wb.active
    
    for row in range(2,ws.max_row+1):
        # Only write out real reactions
        if ws[row] == None or \
                len (ws[row]) == 0 or \
                ws[row][0].value == None:
            continue
        name = ws[row][0].value
        flux_fld = ws[row][flux_col].value
        to_compart = ws[row][to_compart_col].value
        from_compart = ws[row][from_compart_col].value
        dir_fld = ws[row][dir_col].value
        #if not flux_fld is None:
        #    print ("Flux >" + flux_fld + "<")

        compart = from_compart
        if from_compart != to_compart:
            # Its a transport
            compart = from_compart + to_compart
            if compart not in std_names:
                compart = to_compart + from_compart
            if compart not in std_names:
                die ("Illegal compart comb: " + from_compart + to_compart)

        if flux_fld == "max":
            # This is the biomass equation
            if dir_fld == '<':
                name = name + '_rev'
            print (name + "_" + compart, -1.0, 0.0)
        elif flux_fld is None or len(flux_fld) == 0:
            # Unknown flux
            # Output a 0 flux with special 9999 error, indicating unknown
            if dir_fld == '>' or dir_fld == '=':
                print (name + "_" + compart, 0.0, 9999.0)
            # Also/instead put out reverse flow if req'd
            if dir_fld == '<' or dir_fld == '=':
                print (name + '_rev' + "_" + compart, 0.0, 9999.0)
        else:
            fld_tok = flux_fld.split ('±')
            flux = float (fld_tok[0])
            error = float (fld_tok[1])
            if flux < 0:
                # Reversed flow
                flux = -flux
                name = name + '_rev'
            print (name + "_" + compart, flux, error)
                
    
if __name__ == '__main__':
    main()
