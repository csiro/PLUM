#!/usr/bin/env python3
"""TSV Converter for Metabolic Networks.

This module converts TSV (Tab-Separated Values) files containing metabolic
compound and reaction data into a DAT format suitable for flux balance analysis
(FBA) gap-filling tools. It handles compartmentalization, transport reactions,
and reaction reversibility.

The converter reads compound and reaction TSV files, processes stoichiometry
and compartment information, and generates reaction definitions for different
compartment combinations (single, dual, and triple compartment reactions).

Examples
--------
>>> python tsvcvt.py compounds.tsv reactions.tsv output.dat
>>> python tsvcvt.py -o 2 -c ECP compounds.tsv reactions.tsv output.dat

Notes
-----
Compartments are hard-coded by default:
- E: External (input/output layer)
- C: Cytosol
- P: Periplasm
"""

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
    """Metabolic compound representation.

    Represents a metabolic compound with its identifier, name, and usage status.
    Compounds can exist in multiple compartments and are tracked to determine
    which are actually used in reactions.

    Attributes
    ----------
    idstr : str
        Unique identifier for the compound
    name : str
        Human-readable name of the compound
    used : bool
        Flag indicating whether the compound is used in any reaction

    Examples
    --------
    >>> cmpd = Compound('glc', 'Glucose')
    >>> print(cmpd)
    Glucose(glc)
    """
    def __init__(self, idstr, name):
        """Initialize a Compound instance.

        Parameters
        ----------
        idstr : str
            Unique identifier for the compound
        name : str
            Human-readable name of the compound
        """
        self.idstr = idstr
        self.name = name
        self.used = False
        
    def __str__(self):
        """Return string representation of the compound.

        Returns
        -------
        str
            String in format 'name(idstr)'
        """
        return self.name + "(" + self.idstr + ")"

class Reaction:
    """Metabolic reaction representation.

    Represents a metabolic reaction with its identifier, name, reversibility,
    and associated compounds with their stoichiometric coefficients and
    compartment locations.

    Attributes
    ----------
    idstr : str
        Unique identifier for the reaction
    name : str
        Human-readable name of the reaction
    rev : str
        Reversibility indicator
    compound_coeffs : list of list
        List of [compound_id, compartment, coefficient] entries
    max_compart : int
        Maximum compartment index (0=single, 1=dual, 2=triple compartment)

    Examples
    --------
    >>> rxn = Reaction('r001', 'Glucose import', 'reversible')
    >>> rxn.add_compound('glc', 1, -1.0)
    """
    def __init__(self, idstr, name, rev):
        """Initialize a Reaction instance.

        Parameters
        ----------
        idstr : str
            Unique identifier for the reaction
        name : str
            Human-readable name of the reaction
        rev : str
            Reversibility indicator
        """
        self.idstr = idstr
        self.name = name
        self.rev = rev
        self.compound_coeffs = []
        self.max_compart = 0;
        
    def __str__(self):
        """Return string representation of the reaction.

        Returns
        -------
        str
            String in format 'name(idstr)'
        """
        return self.name + "(" + self.idstr + ")"

    def add_compound (self, cmpd, compart, coeff):
        """Add a compound to the reaction.

        Parameters
        ----------
        cmpd : str
            Compound identifier
        compart : int
            Compartment index (0, 1, or 2)
        coeff : float
            Stoichiometric coefficient (negative for reactants, positive for products)

        Notes
        -----
        Updates max_compart to track the highest compartment index used.
        """
        self.compound_coeffs.append ([ cmpd, compart, coeff ])
        if compart > self.max_compart:
            self.max_compart = compart
    
    def is_transport ():
        """Check if the reaction is a transport reaction.

        Returns
        -------
        bool
            True if reaction spans multiple compartments (max_compart > 0)

        Notes
        -----
        This method has a bug: it's missing 'self' parameter and will fail at runtime.
        """
        return self.max_compart > 0

def debug (level, msg):
    """Write debug message to debug file if level threshold is met.

    Parameters
    ----------
    level : int
        Debug level for this message
    msg : str
        Debug message to write

    Notes
    -----
    Messages are only written if level <= debug_level global variable.
    Output is written to the global debug_file.
    """
    if (level <= debug_level):
        print (msg, file=debug_file)


def rmbrace (astr):
    """Remove curly brace expressions from a string.

    Parameters
    ----------
    astr : str
        Input string

    Returns
    -------
    str
        String with all {.*?} patterns removed

    Examples
    --------
    >>> rmbrace('test{remove}string')
    'teststring'
    """
    return re.sub('{.*?}', '', astr)

def rmparen (astr):
    """Remove parentheses from a string.

    Parameters
    ----------
    astr : str
        Input string

    Returns
    -------
    str
        String with all '(' and ')' characters removed

    Examples
    --------
    >>> rmparen('test(remove)string')
    'testremovestring'
    """
    return re.sub('[()]', '', astr)

def rmpbrack (astr):
    """Remove square bracket expressions from a string.

    Parameters
    ----------
    astr : str
        Input string

    Returns
    -------
    str
        String with all [.*?] patterns removed

    Examples
    --------
    >>> rmpbrack('test[remove]string')
    'teststring'
    """
    return re.sub('\[.*?\]', '', astr)

def ends_with (target, str):
    """Check if a string ends with a target substring.

    Parameters
    ----------
    target : str
        Substring to check for at the end
    str : str
        String to check

    Returns
    -------
    bool
        True if str ends with target, False otherwise

    Examples
    --------
    >>> ends_with('ing', 'testing')
    True
    """
    return str[-len(target):] == target

def read_compounds (fn):
    """Read compounds from a TSV file and populate the global compounds dict.

    Parameters
    ----------
    fn : str
        Path to the TSV file containing compound data

    Raises
    ------
    SystemExit
        If expected column headers ('id' in column A, 'name' in column C) are not found

    Notes
    -----
    Expected TSV format:
    - First row: header with 'id' in column 0 and 'name' in column 2
    - Subsequent rows: compound data with id in column 0, name in column 2

    Populates the global 'compounds' dictionary with Compound objects keyed by id.
    """
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
    """Read reactions from a TSV file and populate the global reactions dict.

    Parameters
    ----------
    fn : str
        Path to the TSV file containing reaction data

    Raises
    ------
    SystemExit
        If expected column headers are not found or if a compound referenced
        in stoichiometry is not in the compounds dictionary

    Notes
    -----
    Expected TSV format:
    - Column 0: 'id' - reaction identifier
    - Column 2: 'name' - reaction name
    - Column 4: 'stoichiometry' - format 'coeff:compound_id:compartment' separated by ';'
    - Column 8: 'reversibility'
    - Column 9: 'direction' ('=', '>', or '<')

    For bidirectional reactions ('=' or '<'), creates reversed reaction with '_rev' suffix.
    Populates the global 'reactions' dictionary with Reaction objects.
    """
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
    """Write reactions to a DAT file for flux balance analysis.

    Parameters
    ----------
    fn : str
        Output file path
    compound_fn : str
        Original compound TSV filename (for header comment)
    reaction_fn : str
        Original reaction TSV filename (for header comment)

    Notes
    -----
    Output format:
    - MET lines: Define metabolites for each compound-compartment combination
    - REACTION lines: Define reactions with stoichiometry

    Reactions are expanded based on compartment configuration:
    - max_compart = 0: Single compartment reactions (E, C, or P)
    - max_compart = 1: Dual compartment transport (EP, PC, CE)
    - max_compart = 2: Triple compartment transport (EPC, PCE)

    Uses global variables: compartments, default_obj_coeff, bad_obj_coeff
    Reads from global dictionaries: compounds, reactions
    """
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
    """Print error message to stderr and exit with status 1.

    Parameters
    ----------
    message : str
        Error message to display

    Raises
    ------
    SystemExit
        Always exits with code 1
    """
    print(message, file=sys.stderr)
    exit (1)

def usage (str = ""):
    """Print usage information and exit.

    Parameters
    ----------
    str : str, optional
        Additional message to print before usage information (default: '')

    Raises
    ------
    SystemExit
        Always exits with code 1

    Notes
    -----
    Displays command syntax, argument descriptions, and available switches.
    """
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
    """Main entry point for TSV to DAT conversion.

    Parses command-line arguments, reads compound and reaction TSV files,
    and generates a DAT file for flux balance analysis.

    Command-line Arguments
    ----------------------
    compounds.tsv : str
        Path to TSV file containing compound definitions
    reactions.tsv : str
        Path to TSV file containing reaction definitions
    out_fn : str, optional
        Output DAT file path (default: 'out.dat')

    Options
    -------
    -o # : int
        Objective coefficient for reactions (default: 1)
    -c str : str
        Compartment letters, one per compartment (default: 'CPE')
    -a str : str
        Default compartments for reactions (default: 'CP')

    Raises
    ------
    SystemExit
        If insufficient arguments provided or parsing errors occur

    Notes
    -----
    Modifies global variables: compartments, default_compart_names, default_obj_coeff
    Populates global dictionaries: compounds, reactions
    """
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
