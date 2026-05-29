#!/usr/bin/env python3

import sys
import os
import re
import fileinput
import xml.dom.minidom as MD
import xml.etree.ElementTree as ET
from datetime import datetime


debug_level = 3
debug_file = open ("debug.out", "w")

default_obj_coeff = 2000
gene_ind_obj_coeff = 1
include_ex = True
include_dm = True
allow_dm_rev = False
out_fn = "out.pld"
fix_name = True

gene_products = {}
species = {}
reactions = {}

obj_coeffs = {}

param = {}

skip_ex_count = 0
skip_dm_count = 0
fix_dm_count = 0

class GeneProduct:
    
    def __init__(self):
        self.name = ''
        self.idstr = ''
        
    def __str__(self):
        return self.name + "(" + self.idstr + ")"
    
class Species:
    def __init__(self):
        self.name = ''
        self.idstr = ''
        
    def __str__(self):
        return self.name + "(" + self.idstr + ")"
    
class Reaction:
    def __init__(self):
        global default_obj_coeff
        self.name = ''
        self.idstr = ''
        self.reversible = False
        self.obj_coeff = default_obj_coeff
        self.lowerFluxBound = 0
        self.upperFluxBound = 0
        self.species_coeffs = []
        self.gene_indicated = False
        
    def __str__(self):
        return self.name + "(" + self.idstr + ")"
    
    def add_species (self, spec, coeff):
        self.species_coeffs.append ([ spec, coeff ])
        
    def is_ex(self):
        if len(self.species_coeffs) == 0:
            return False
        # EXs have no inputs
        for sp, coeff in self.species_coeffs:
            if coeff < 0:
                return False
        return True
        
    def is_dm(self):
        if len(self.species_coeffs) == 0:
            return False
        # DMs have no outputs
        for sp, coeff in self.species_coeffs:
            if coeff > 0:
                return False
        return True
        
    def set_gene_indicated (self, val):
        self.gene_indicated = val
        if val:
            self.obj_coeff = gene_ind_obj_coeff
        else:
            self.obj_coeff = default_obj_coeff
    
    def coeff_str (self):
        retstr = ''
        sep = ''
        for sp, coeff in self.species_coeffs:
            if coeff < 0:
                retstr += sep
                if coeff != -1:
                    retstr += str(-coeff) + ' '
                retstr += sp.idstr
                sep = ' + '
        retstr += ' -> '
        sep = ' '
        for sp, coeff in self.species_coeffs:
            if coeff > 0:
                retstr += sep
                if coeff != 1:
                    retstr += str(coeff) + ' '
                retstr += sp.idstr
                sep = ' + '
        return retstr


def debug (level, msg):
    if (level <= debug_level):
        print (msg, file=debug_file)


def rmbrace (astr):
    return re.sub('{.*?}', '', astr)

def ends_with (target, str):
    return str[-len(target):] == target

def proc_root (myroot):
    for elt in myroot[0]:
        tag = rmbrace(elt.tag)
        debug (2, "Root tag " + tag)
        if tag == 'listOfGeneProducts':
            proc_listOfGeneProducts (elt)
        if tag == 'listOfObjectives':
            proc_listOfObjectives (elt)
        if tag == 'listOfSpecies':
            proc_listOfSpecies (elt)
        if tag == 'listOfParameters':
            proc_listOfParameters (elt)
        if tag == 'listOfReactions':
            proc_listOfReactions (elt)

def proc_listOfGeneProducts (head):
    debug (2, "proc_listOfGeneProducts")
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "  Tag " + tag)
        if tag == 'geneProduct':
            proc_geneProduct (elt)
    
def proc_listOfObjectives (head):
    debug (2, "proc_listOfObjectives")
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "  Tag " + tag)
        if tag == 'objective':
            proc_objective (elt)
    
def proc_objective (head):
    debug (2, "proc_objective")
    idstr = ''
    for a in head.attrib:
        attrib = rmbrace (a)
        val = rmbrace (head.attrib[a])
        debug (3, "    Attrib " + attrib + " value " + val)
        if attrib == "id":
            idstr = val
        elif attrib == "type":
            if val != 'maximize':
                print ("Can only handle maximisation objectives - objective " + idstr)
                debug_file.close()
                exit(2)
                
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (2, "  Tag " + tag)
        if tag == 'listOfFluxObjectives':
            proc_listOfFluxObjectives (elt)
    
def proc_listOfFluxObjectives (head):
    debug (2, "proc_listOfFluxObjectives")
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (2, "  Tag " + tag)
        if tag == 'fluxObjective':
            proc_fluxObjective (elt)
        else:
            print ("Can only handle flux objectives")
            debug_file.close()
            exit(2)
        
    
def proc_fluxObjective (head):
    debug (2, "proc_fluxObjectives")
    idstr = ''
    for a in head.attrib:
        attrib = rmbrace (a)
        val = rmbrace (head.attrib[a])
        debug (3, "    Attrib " + attrib + " value " + val)
        if attrib == "reaction":
            idstr = val
        elif attrib == "coefficient":
            obj_coeffs[idstr] = val
            print ("Objective reaction is " + idstr)
            
def proc_listOfSpecies (head):
    debug (2, "proc_listOfSpecies")
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "  Tag " + tag)
        if tag == 'species':
            proc_species (elt)
    
def proc_listOfParameters (head):
    debug (2, "proc_listOfParameters")
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "  Tag " + tag)
        if tag == 'parameter':
            proc_parameter (elt)
            
def proc_listOfReactions (head):
    debug (2, "proc_listOfReactions")
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "  Tag " + tag)
        if tag == 'reaction':
            proc_reaction (elt)
    
def proc_geneProduct (head):
    gp = GeneProduct()
    for a in head.attrib:
        attrib = rmbrace (a)
        val = rmbrace (head.attrib[a])
        debug (3, "    Attrib " + attrib + " value " + val)
        if attrib == "name":
            gp.name = val
        elif attrib == "id":
            gp.idstr = val

    if gp.idstr in gene_products:
        print ('Duplicate gene product: ' + idstr)
        debug_file.close()
        exit(4)
        
    gene_products[gp.idstr] = gp

def proc_species (head):
    s = Species()
    name = ""
    formula = ""
    for a in head.attrib:
        attrib = rmbrace (a)
        val = rmbrace (head.attrib[a])
        debug (3, "    Attrib " + attrib + " value " + val)
        if attrib == "name":
            name = val
        elif attrib == "id":
            s.idstr = val
        elif attrib == "chemicalFormula":
            formula = val
            
    if formula == "":
        s.name = name
    else:
        s.name = name + " [" + formula + "]"
            
    if s.idstr in species:
        print ('Duplicate species: ' + idstr)
        debug_file.close()
        exit(4)

    species[s.idstr] = s
        
def proc_parameter (head):
    name = ''
    value = ''
    for a in head.attrib:
        attrib = rmbrace (a)
        val = rmbrace (head.attrib[a])
        debug (3, "    Attrib " + attrib + " value " + val)
        if attrib == "id":
            name = val
        elif attrib == "value":
            param[name] = val
            
def proc_reaction (head):
    global include_ex
    global include_dm
    global skip_ex_count
    global skip_dm_count
    global fix_dm_count
    global allow_dm_rev
    
    react = Reaction()
    for a in head.attrib:
        attrib = rmbrace (a)
        val = rmbrace (head.attrib[a])
        debug (3, "    Attrib " + attrib + " value " + val)
        if attrib == "name":
            react.name = val
        elif attrib == "id":
            react.idstr = val
        elif attrib == "lowerFluxBound":
            if val in param:
                react.lowerFluxBound = float(param[val])
            else:
                react.lowerFluxBound = float(val)
        elif attrib == "upperFluxBound":
            if val in param:
                react.upperFluxBound = float(param[val])
            else:
                react.upperFluxBound = float(val)
        elif attrib == "reversible":
            react.reversible = (val == "true")
            
    # Also process sub-elements            
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "    Tag " + tag)
        if tag == 'listOfReactants':
            proc_listOfReactants (elt, react)
        elif tag == 'listOfProducts':
            proc_listOfProducts (elt, react)
        elif tag == 'geneProductAssociation':
            proc_geneProductAssociation (elt, react)
            
    if react.idstr in reactions.keys():
        print ('Duplicate reaction: ' + idstr)
        debug_file.close()
        exit(4)
          
    debug (2, "  React " + react.name + " len " + str(len(react.species_coeffs)))
    if react.is_ex() and not include_ex:
        print ('Skip ex reaction ' + react.name)
        skip_ex_count += 1
        return
    if react.is_dm() and not include_dm:
        print ('Skip dm reaction ' + react.name)
        skip_dm_count += 1
        return
    
    if not allow_dm_rev and react.is_dm() and \
            react.reversible and react.lowerFluxBound < 0:
        debug ( \
            2, 'Make DM react ' + react.idstr + ' ' + react.name + \
            ' non-reversible')
        debug (2, '  ' + react.coeff_str())
        fix_dm_count += 1
        react.lowerFluxBound = 0
        
    reactions[react.idstr] = react
        
def proc_listOfReactants (head, react):
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "      Tag " + tag)
        if tag == 'speciesReference':
            proc_speciesReference(elt, react, -1.0)

def proc_listOfProducts (head, react):
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "      Tag " + tag)
        if tag == 'speciesReference':
            proc_speciesReference(elt, react, +1.0)

def proc_speciesReference (head, react, sign):
    idstr = ''
    coeff = ''
    for a in head.attrib:
        attrib = rmbrace (a)
        val = rmbrace (head.attrib[a])
        debug (3, "          Attrib " + attrib + " value " + val)
        if attrib == "species":
            idstr = val
        elif attrib == "stoichiometry":
            coeff = float(val) * sign
            
    if not idstr in species:
        print ("Lost species id " + idstr + " in reaction " + react.idstr)
        debug_file.close()
        exit(1)

    react.add_species (species[idstr], coeff)
    debug (
        2, "          React " + react.idstr + " coeff count " + 
        str(len(react.species_coeffs)))

            
def proc_geneProductAssociation (head, react):
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "      Tag " + tag)
        if tag == 'or':
            proc_geneProductAssociation (elt, react)
        elif tag == 'and':
            proc_geneProductAssociation (elt, react)
        elif tag == 'geneProductRef':
            react.set_gene_indicated (True)

def write_pld_file (fn, xml_fn):
    global include_ex
    global include_dm
    global allow_dm_rev
    
    dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S") 

    with open (fn, "w") as out:
        print ('# Created ' + dt_string + ' by sbmlcvt.py from ' + xml_fn, file=out)
        print (
            '# Using default_obj_coeff ' + str(default_obj_coeff) + 
            ' gene_ind_obj_coeff ' + str(gene_ind_obj_coeff) + 
            ' include_ex ' + str(include_ex) +
            ' include_dm ' + str(include_dm) +
            ' allow_dm_rev ' + str(allow_dm_rev), file=out)
        
        for spec in species.values():
            if fix_name:
                spec.idstr = re.sub ("^M_", "", spec.idstr)
            print ('MET ' + spec.idstr + '	"' + spec.name + '"', file=out)
    
        for react in reactions.values():
            if fix_name:
                react.idstr = re.sub ("^R_", "", react.idstr)
            if react.upperFluxBound > 0:
                # Remove R_ if needed
                print ( \
                    'REACTION ' + react.idstr, \
                        react.obj_coeff, \
                        react.upperFluxBound, \
                        file=out, end=' ')
                for sp, coeff in react.species_coeffs:
                    print (sp.idstr + ' ' + str(coeff), file=out, end=' ')
                if len(react.name) > 0:
                    print (" __fullname__ \"" + react.name, file=out, end='\"')
                print (file=out)

            if react.reversible:
                # Output with negated ceoffs
                # Use negated lower bound (which should be -ve) as UB
                if react.lowerFluxBound < 0:
                    print ( \
                        'REACTION ' + react.idstr + "__rev", \
                            react.obj_coeff, \
                            -react.lowerFluxBound, \
                            file=out, end=' ')
                    for sp, coeff in react.species_coeffs:
                        print (sp.idstr + ' ' + str(-coeff), file=out, end=' ')
                    if len(react.name) > 0:
                        print (" __fullname__ \"" + react.name + " (reverse)", file=out, end='\"')
                    print (file=out)
    
    print ("Done. wrote output to " + fn)
    
def die (message):
    print(message, file=sys.stderr)
    exit (1)

def usage (msg = ""):
    print (msg)
    print (sys.argv[0], end=' ')
    print ('''[-j #] [-g #] [-j #] [-ex] [-dm] [-exlb] [+dmrev] [-rm] file.xml [out_fn] 
        Parse and convert an SBML file on metabolic reactions
    Args
        file.xml: file to process
	      out_fn: Output file [''' + out_fn + ''']
    Switches
              -j: default reaction objective coeff [''' + str(default_obj_coeff) + ''']
              -g: gene-indicated reaction objective coeff [''' + str(gene_ind_obj_coeff) + ''']
             -rm: don't remove R_ and M_ from react/met names [do]
             -ex: don't include exchange reactions [do]
                  (exchange reactions are those that have 
                   a single product, and no reactants)
             -dm: don't include demand reactions [do]
                  (demand reactions are those that have 
                   a single reactant, and no products)
          +dmrev: Allow reversible DM reactions [don't]
                  By default, DM reactions (one reactant, no products) 
                  will NOT be allowed to reverse. This switch means
                  they will be allowed to supply the metabolite.
''')
    exit (1)


def main():
    # parse argv
    global default_obj_coeff
    global gene_ind_obj_coeff
    global include_ex
    global include_dm
    global allow_dm_rev

    p = 0
    xml_fn = ""
    out_fn = ""
    
    while len(sys.argv) > 1:
        arg = sys.argv.pop(1)
        if arg[0] == '-' or arg[0] == '+':
            if arg == '-j' and len(sys.argv) > 1:
                default_obj_coeff = float(sys.argv.pop(1))
            if arg == '-g' and len(sys.argv) > 1:
                gene_ind_obj_coeff = float(sys.argv.pop(1))
            elif arg == '-ex':
                include_ex = False
            elif arg == '-dm':
                include_dm = False
            elif arg == '+dmrev':
                allow_dm_rev = True
            elif arg == '-rm':
                fix_name = False
            else:
                usage("Unrecognised arg: " + arg)
        else:
            # Alternatively
            if xml_fn == "":
                xml_fn = arg
            elif out_fn == "":
                out_fn = arg
            else:
                usage("Too many file names: " + arg + " xml " + xml_fn + " out " + out_fn)
                
    # Or, do stuff with named file
    if xml_fn == "":
        usage("Need filename")
        
    if out_fn == "":
        out_fn = "out.pld"

    xmldoc = ET.parse (xml_fn)
    
    myroot = xmldoc.getroot()
    
    proc_root (myroot)
    
    for react, val in obj_coeffs.items():
        if not react in reactions:
            print ("Lost objective reaction", react)
            exit(3)
        reactions[react].obj_coeff = -1.0
        print ("Found objective reaction " + react)
        
    debug (1, "List of params")
    for key, val in param.items():
        debug (1, "  " + key + ": " + val)
        
    debug (1, "List of gene products:")
    for gp in gene_products.values():
        debug (1, "  " + str(gp))
    debug (1, "List of species:")   
    for s in species.values():
        debug (1, "  " + str(s))

    debug (1, "List of reactions:")   
    for react in reactions.values():
        debug (1, "  " + react.idstr)
        debug (1, "     name " + react.name)
        debug (1 , "      obj " + str(react.obj_coeff))
        debug (1, "       lb " + str(react.lowerFluxBound))
        debug (1, "       ub " + str(react.upperFluxBound))
        debug (1, "      rev " + str(react.reversible))
        debug (1, "     species coeffs ")
        for sp, coeff in react.species_coeffs:
            debug (1, "          " + sp.idstr + ': ' + str(coeff))
            
    write_pld_file (out_fn, xml_fn)
    
    if not include_ex:
        print ('Skipped', skip_ex_count, 'EX reactions')
    if not include_dm:
        print ('Skipped', skip_dm_count, 'DM/sink reactions')
    if not allow_dm_rev:
        print ('Updated', fix_dm_count, 'lower bounds on DM/sink reactions')
    
        
    print('All done')

if __name__ == '__main__':
    main()
