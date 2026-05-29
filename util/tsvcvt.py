#!/usr/bin/env python3

import sys
import os
import re
import fileinput
from datetime import date


debug_level = 6
debug_file = open ("debug.out", "w")

default_obj_coeff = 1
bad_obj_coeff = 9999
# Compartments are hard-coded
# "External" (E) is our input/output layer, typically moving through 
# Periplasm (P) to # Cytosol (C), and then to E again
compartments = "ECP" # External, Cytosol, Periplasm

# Compartments for non-transport reactions
react1_comparts = { "E", "C", "P" } 
# Compartments for transport reactions one compart to another
react2_comparts = { "EP", "PC", "CE" } 
# Compartments for transport reactions across 3 compartments
react3_comparts = { "EPC", "PCE"}

compounds = {}
reactions = {}

class Compound:
    def __init__(self, idstr, name):
        self.idstr = idstr
        self.name = name
        self.used = False
        
    def __str__(self):
        return self.name + "(" + self.idstr + ")"
    
class Reaction:
    def __init__(self, idstr, name, rev):
        self.idstr = idstr
        self.name = name
        self.rev = rev
        self.compound_coeffs = []
        self.max_compart = 0;
        
    def __str__(self):
        return self.name + "(" + self.idstr + ")"
    
    def add_compound (self, cmpd, compart, coeff):
        self.compound_coeffs.append ([ cmpd, compart, coeff ])
        if compart > self.max_compart:
            self.max_compart = compart
    
    def is_transport ():
        return self.max_compart > 0

def debug (level, msg):
    if (level <= debug_level):
        print (msg, file=debug_file)


def rmbrace (astr):
    return re.sub('{.*?}', '', astr)

def rmparen (astr):
    return re.sub('[()]', '', astr)

def rmpbrack (astr):
    return re.sub('\[.*?\]', '', astr)

def ends_with (target, str):
    return str[-len(target):] == target

def read_compounds (fn):
    count = 0
    with open(fn) as infile:
        for line in infile:
            count += 1
            fld = line.rstrip().split("\t")
            if count == 1:
                # Check field names
                if fld[0] != 'id':
                    die ("Col A expected 'id', got " + fld[0])
                if fld[2] != 'name':
                    die ("Col C expected 'name', got " + fld[2])
            else:
                idstr = fld[0]
                name = fld[2]
                cmpd = Compound (idstr, name)
                compounds[idstr] = cmpd

def read_reactions (fn):
    count = 0
    with open(fn) as infile:
        for line in infile:
            count += 1
            fld = line.rstrip().split("\t")
            if count == 1:
                # Check field names
                if fld[0] != 'id':
                    die ("Col A expected 'id', got " + fld[0])
                if fld[2] != 'name':
                    die ("Col C expected 'name', got " + fld[2])
                if fld[4] != 'stoichiometry':
                    die ("Col G expected 'stoichiometry', got " + fld[4])
                if fld[8] != 'reversibility':
                    die ("Col I expected 'reversibility', got " + fld[8])
                if fld[9] != 'direction':
                    die ("Col J expected 'direction', got " + fld[9])
            else:
                idstr = fld[0]
                name = fld[2]
                stoich = fld[4]
                rev = fld[8]
                direction = fld[9]
                react = Reaction (idstr, name, rev)
                if direction == '=' or direction == '>':
                    reactions[idstr] = react
                    
                fld = stoich.split(";")
                debug (5, "Line " +  str(count) + " Stoich: " + stoich)
                print (idstr + " " + stoich)
                for colsepstr in fld:
                    debug (5, "  colsepstr: " + colsepstr)
                    tok = colsepstr.split(":")
                    if len(tok) > 1:
                        coeff = float(tok[0].replace('"',''))
                        cmpd_id = tok[1]
                        compart = int(tok[2])
                            
                        debug (6, "      coeff: " + str(coeff))
                        debug (6, "       cmpd: " + cmpd_id)
                        debug (6, "    compart: " + str(compart))
                        if not cmpd_id in compounds:
                            die ("Lost compound id " + cmpd_id + \
                                     " in reaction " + react.idstr)
                        cmpd = compounds[cmpd_id]
                        cmpd.used = True
                            
                        react.add_compound (cmpd_id, compart, coeff)
                        
                if direction == '=' or direction == '<':
                    # Add the reverse reaction
                    idstr = idstr + '_rev'
                    name = name + ' reversed'
                    rev_react = \
                        Reaction (idstr, name, rev)
                    reactions[idstr] = rev_react
                    for cmpd_id, compart, coeff in react.compound_coeffs:
                        rev_react.add_compound (cmpd_id, compart, -coeff)
                        
            
def write_datfile (fn, compound_fn, reaction_fn):  
    global compartments
    global default_obj_coeff
    global bad_obj_coeff
    
    dt_string = date.today().strftime("%d/%m/%Y %H:%M:%S") 

    with open (fn, "w") as out:
        print (\
            '# Created by tsvcvt.py on ' + dt_string + \
                " from " + compound_fn + ", "  + reaction_fn, \
                file=out)

        for cmpd in compounds.values():
            if cmpd.used:
                for c in compartments:
                    print ('MET ' + cmpd.idstr + "_" + c, file=out, end='')
                    print ('	# ' + cmpd.name, file=out)
                    
        for react in reactions.values():
            if react.max_compart == 0:
                # Not a transport 
                for c in react1_comparts:
                    print ( \
                        'REACTION ' + react.idstr + "_" + c + " " +  \
                            str(default_obj_coeff), \
                            file=out, end=' ')
                    for sp, compart, coeff in react.compound_coeffs:
                        print (sp + "_" + c + " ", coeff, file=out, end=' ')
                    print ('  # ' + react.name + " " + c, file=out)
            elif react.max_compart == 1:
                for str2 in react2_comparts:
                    print (
                        'REACTION ' + react.idstr + "_" + str2 + " " + 
                        str(default_obj_coeff), 
                        file=out, end=' ')
                    for sp, compart, coeff in react.compound_coeffs:
                        print ( \
                            sp + "_" + str2[compart] + " ", \
                                coeff, file=out, end=' ')
                    print (\
                        '  # ' + react.name + " transp " + str2, \
                            file=out)
            elif react.max_compart == 2:
                for str3 in react3_comparts:
                    print ( \
                        'REACTION ' + react.idstr + "_" + str3 + " " + \
                            str(default_obj_coeff), \
                            file=out, end=' ')
                    for sp, compart, coeff in react.compound_coeffs:
                        print ( \
                            sp + "_" + str3[compart] + " ", coeff, \
                                file=out, end=' ')
                    print ( \
                        '  # ' + react.name + " transp " + str3, \
                            file=out)


def die (message):
    print(message, file=sys.stderr)
    exit (1)

def usage (str = ""):
    print (str)
    print (sys.argv[0], end=' ')
    print ('''[-o #] [-c str] compounds.tsv reactions.tsv [out_fn] 
        Parse and convert an TSV file on metabolic reactions
    Args
        compounds.tsv: Compounds
         reaction.tsv: Reactions
	           out_fn: Output file [out.dat]
    Switches
        -o: objective coefficient [1]
        -c: Compartment letters. One letter for each compartment [CPE]
        -a: Default compartments for reactions [CP]
''')
    exit (1)


def main():
    global compartments
    global default_compart_names
    global default_obj_coeff
    
    # parse argv
    compound_fn = ""
    reaction_fn = ""
    out_fn = ""
    
    while len(sys.argv) > 1:
        arg = sys.argv.pop(1)
        if arg[0] == '-':
            if arg == '-o' and len(sys.argv) > 1:
                default_obj_coeff = int(sys.argv.pop(1))
            elif arg == '-c' and len(sys.argv) > 1:
                compartments = sys.argv.pop(1)
            elif arg == '-a' and len(sys.argv) > 1:
                default_compart_names = sys.argv.pop(1)
            else:
                usage("Unrecognised arg: " + arg)
        else:
            # Alternatively
            if compound_fn == "":
                compound_fn = arg
            elif reaction_fn == "":
                reaction_fn = arg
            elif out_fn == "":
                out_fn = arg
            else:
                usage("Too many file names: " + arg)
                
    if reaction_fn == "":
        usage("Need at least 2 filenames")
    if out_fn == "":
        out_fn = "out.dat"
        
    read_compounds (compound_fn);    
    read_reactions (reaction_fn);

    debug (1, "List of compounds:")   
    for s in compounds.values():
        debug (1, "  " + str(s))

    debug (1, "List of reactions:")   
    for react in reactions.values():
        debug (1, "  " + react.idstr)
        debug (1, "     name " + react.name)
        debug (1, "       rev " + react.rev)
        debug (1, "     compound coeffs ")
        for sp, compart, coeff in react.compound_coeffs:
            debug (1, "          " + sp + ': ' + str(coeff))
            
    write_datfile (out_fn, compound_fn, reaction_fn)
        
    print('Wrote ' + out_fn)

if __name__ == '__main__':
    main()
