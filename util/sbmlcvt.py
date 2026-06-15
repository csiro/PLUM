#!/usr/bin/env python3
"""SBML to PLD converter for metabolic flux balance analysis.

This module parses Systems Biology Markup Language (SBML) files containing
metabolic reaction networks and converts them to PLD (PLUM data) format for
use in gap-filling and flux balance analysis tools. It handles gene products,
species (metabolites), reactions with stoichiometry, flux bounds, and objective
coefficients.

The converter supports filtering of exchange (EX) and demand (DM) reactions,
customization of objective coefficients, and automatic handling of reversible
reactions by creating forward and reverse reaction entries.

Global Configuration
--------------------
debug_level : int
    Verbosity level for debug output (default: 3)
default_obj_coeff : int
    Default objective coefficient for reactions (default: 2000)
gene_ind_obj_coeff : int
    Objective coefficient for gene-indicated reactions (default: 1)
include_ex : bool
    Whether to include exchange reactions (default: True)
include_dm : bool
    Whether to include demand reactions (default: True)
allow_dm_rev : bool
    Whether to allow reversible demand reactions (default: False)
fix_name : bool
    Whether to strip R_ and M_ prefixes from names (default: True)

Examples
--------
Convert an SBML file with default settings:

    $ python sbmlcvt.py model.xml output.pld

Exclude exchange reactions and set custom objective coefficient:

    $ python sbmlcvt.py -j 1000 -ex model.xml output.pld

Notes
-----
The output PLD format includes:
- MET entries: metabolite ID and name with chemical formula
- REACTION entries: reaction ID, objective coefficient, upper flux bound,
  stoichiometry, and optional full name
- Reversible reactions are split into forward and reverse entries
"""

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
    """Represents a gene product in the metabolic model.

    Gene products are associated with reactions through gene product associations,
    indicating genetic evidence for the reaction's occurrence in the organism.

    Attributes
    ----------
    name : str
        Human-readable name of the gene product
    idstr : str
        Unique identifier string for the gene product

    Examples
    --------
    >>> gp = GeneProduct()
    >>> gp.name = "b0001"
    >>> gp.idstr = "G_b0001"
    >>> str(gp)
    'b0001(G_b0001)'
    """

    def __init__(self):
        """Initialize a GeneProduct instance.

        Creates a new gene product with empty name and identifier strings.
        """
        self.name = ''
        self.idstr = ''
        
    def __str__(self):
        """Return string representation of the gene product.

        Returns
        -------
        str
            Formatted string as 'name(idstr)'
        """
        return self.name + "(" + self.idstr + ")"

class Species:
    """Represents a metabolic species (metabolite) in the model.

    Species participate in reactions as reactants or products with associated
    stoichiometric coefficients.

    Attributes
    ----------
    name : str
        Human-readable name, optionally including chemical formula
    idstr : str
        Unique identifier string for the species

    Examples
    --------
    >>> species = Species()
    >>> species.name = "Glucose [C6H12O6]"
    >>> species.idstr = "M_glc__D_c"
    >>> str(species)
    'Glucose [C6H12O6](M_glc__D_c)'
    """
    def __init__(self):
        """Initialize a Species instance.

        Creates a new species with empty name and identifier strings.
        """
        self.name = ''
        self.idstr = ''

    def __str__(self):
        """Return string representation of the species.

        Returns
        -------
        str
            Formatted string as 'name(idstr)'
        """
        return self.name + "(" + self.idstr + ")"

class Reaction:
    """Represents a metabolic reaction in the flux balance analysis model.

    Reactions define transformations of metabolites with stoichiometric coefficients,
    flux bounds, reversibility, and objective function coefficients. Reactions may be
    classified as exchange (EX) or demand (DM) reactions based on their stoichiometry.

    Attributes
    ----------
    name : str
        Human-readable name of the reaction
    idstr : str
        Unique identifier string for the reaction
    reversible : bool
        Whether the reaction can proceed in reverse direction
    obj_coeff : int or float
        Objective function coefficient for optimization
    lowerFluxBound : float
        Lower bound on reaction flux (typically negative for reversible reactions)
    upperFluxBound : float
        Upper bound on reaction flux
    species_coeffs : list of [Species, float]
        List of (species, coefficient) pairs; negative for reactants, positive for products
    gene_indicated : bool
        Whether the reaction has gene product association evidence

    Examples
    --------
    >>> rxn = Reaction()
    >>> rxn.name = "Glucose transport"
    >>> rxn.idstr = "R_GLCt"
    >>> rxn.reversible = False
    >>> rxn.add_species(glucose, -1.0)
    >>> rxn.add_species(glucose_internal, 1.0)
    """
    def __init__(self):
        """Initialize a Reaction instance.

        Creates a new reaction with default values including the global default
        objective coefficient and empty species list.
        """
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
        """Return string representation of the reaction.

        Returns
        -------
        str
            Formatted string as 'name(idstr)'
        """
        return self.name + "(" + self.idstr + ")"

    def add_species (self, spec, coeff):
        """Add a species with its stoichiometric coefficient to the reaction.

        Parameters
        ----------
        spec : Species
            The metabolite species participating in the reaction
        coeff : float
            Stoichiometric coefficient (negative for reactants, positive for products)
        """
        self.species_coeffs.append ([ spec, coeff ])
        
    def is_ex(self):
        """Check if this is an exchange (EX) reaction.

        Exchange reactions have only products (no reactants) and typically represent
        metabolite import from the environment.

        Returns
        -------
        bool
            True if the reaction has no negative coefficients (no inputs), False otherwise
        """
        if len(self.species_coeffs) == 0:
            return False
        # EXs have no inputs
        for sp, coeff in self.species_coeffs:
            if coeff < 0:
                return False
        return True
        
    def is_dm(self):
        """Check if this is a demand (DM) reaction.

        Demand reactions have only reactants (no products) and typically represent
        metabolite consumption or export to biomass/maintenance.

        Returns
        -------
        bool
            True if the reaction has no positive coefficients (no outputs), False otherwise
        """
        if len(self.species_coeffs) == 0:
            return False
        # DMs have no outputs
        for sp, coeff in self.species_coeffs:
            if coeff > 0:
                return False
        return True
        
    def set_gene_indicated (self, val):
        """Set gene indication status and update objective coefficient accordingly.

        Gene-indicated reactions have genetic evidence and receive a different
        (typically lower) objective coefficient to prioritize them in gap-filling.

        Parameters
        ----------
        val : bool
            True if the reaction has gene product association, False otherwise
        """
        self.gene_indicated = val
        if val:
            self.obj_coeff = gene_ind_obj_coeff
        else:
            self.obj_coeff = default_obj_coeff
    
    def coeff_str (self):
        """Generate a human-readable stoichiometry string for the reaction.

        Creates a string representation in the format:
        'coeff1 reactant1 + coeff2 reactant2 -> coeff3 product1 + coeff4 product2'
        Coefficients of 1 or -1 are omitted.

        Returns
        -------
        str
            Formatted stoichiometry string with reactants, arrow, and products
        """
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
    """Write debug message to debug file if level threshold is met.

    Parameters
    ----------
    level : int
        Debug priority level of this message
    msg : str
        Debug message to write

    Notes
    -----
    Messages are only written if level <= debug_level (global variable).
    Output goes to the global debug_file handle.
    """
    if (level <= debug_level):
        print (msg, file=debug_file)


def rmbrace (astr):
    """Remove XML namespace declarations in braces from a tag string.

    Parameters
    ----------
    astr : str
        XML tag or attribute string that may contain namespace in {namespace} format

    Returns
    -------
    str
        String with all {.*} patterns removed

    Examples
    --------
    >>> rmbrace('{http://www.sbml.org/sbml/level3/version1/core}reaction')
    'reaction'
    """
    return re.sub('{.*?}', '', astr)

def ends_with (target, str):
    """Check if a string ends with the specified target suffix.

    Parameters
    ----------
    target : str
        Suffix string to check for
    str : str
        String to test

    Returns
    -------
    bool
        True if str ends with target, False otherwise
    """
    return str[-len(target):] == target

def proc_root (myroot):
    """Process the root element of the SBML XML tree.

    Dispatches to appropriate list processors based on XML tag names.
    Handles listOfGeneProducts, listOfObjectives, listOfSpecies, listOfParameters,
    and listOfReactions.

    Parameters
    ----------
    myroot : xml.etree.ElementTree.Element
        Root element of the parsed SBML XML document

    Notes
    -----
    Modifies global dictionaries: gene_products, species, reactions, param
    """
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
    """Process the listOfGeneProducts element from SBML.

    Iterates through all geneProduct child elements and processes each one.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The listOfGeneProducts XML element

    Notes
    -----
    Populates the global gene_products dictionary.
    """
    debug (2, "proc_listOfGeneProducts")
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "  Tag " + tag)
        if tag == 'geneProduct':
            proc_geneProduct (elt)
    
def proc_listOfObjectives (head):
    """Process the listOfObjectives element from SBML.

    Iterates through all objective child elements and processes each one.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The listOfObjectives XML element

    Notes
    -----
    Populates the global obj_coeffs dictionary with objective reaction coefficients.
    """
    debug (2, "proc_listOfObjectives")
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "  Tag " + tag)
        if tag == 'objective':
            proc_objective (elt)
    
def proc_objective (head):
    """Process a single objective element from SBML.

    Extracts objective attributes and processes flux objectives. Only maximization
    objectives are supported.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The objective XML element

    Raises
    ------
    SystemExit
        If the objective type is not 'maximize'

    Notes
    -----
    Delegates to proc_listOfFluxObjectives for child elements.
    """
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
    """Process the listOfFluxObjectives element from an SBML objective.

    Iterates through fluxObjective child elements. Only flux objectives are supported.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The listOfFluxObjectives XML element

    Raises
    ------
    SystemExit
        If any non-fluxObjective child element is encountered
    """
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
    """Process a single fluxObjective element from SBML.

    Extracts the reaction ID and objective coefficient for the flux objective.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The fluxObjective XML element

    Notes
    -----
    Populates the global obj_coeffs dictionary with reaction ID as key.
    """
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
    """Process the listOfSpecies element from SBML.

    Iterates through all species child elements and processes each one.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The listOfSpecies XML element

    Notes
    -----
    Populates the global species dictionary.
    """
    debug (2, "proc_listOfSpecies")
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "  Tag " + tag)
        if tag == 'species':
            proc_species (elt)
    
def proc_listOfParameters (head):
    """Process the listOfParameters element from SBML.

    Iterates through all parameter child elements and processes each one.
    Parameters typically define flux bounds and other model constants.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The listOfParameters XML element

    Notes
    -----
    Populates the global param dictionary.
    """
    debug (2, "proc_listOfParameters")
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "  Tag " + tag)
        if tag == 'parameter':
            proc_parameter (elt)
            
def proc_listOfReactions (head):
    """Process the listOfReactions element from SBML.

    Iterates through all reaction child elements and processes each one.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The listOfReactions XML element

    Notes
    -----
    Populates the global reactions dictionary.
    """
    debug (2, "proc_listOfReactions")
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "  Tag " + tag)
        if tag == 'reaction':
            proc_reaction (elt)
    
def proc_geneProduct (head):
    """Process a single geneProduct element from SBML.

    Extracts gene product name and ID attributes and stores in global dictionary.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The geneProduct XML element

    Raises
    ------
    SystemExit
        If a duplicate gene product ID is encountered

    Notes
    -----
    Adds the gene product to the global gene_products dictionary.
    """
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
    """Process a single species element from SBML.

    Extracts species name, ID, and chemical formula. If a formula is present,
    it is appended to the name in brackets.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The species XML element

    Raises
    ------
    SystemExit
        If a duplicate species ID is encountered

    Notes
    -----
    Adds the species to the global species dictionary.
    """
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
    """Process a single parameter element from SBML.

    Extracts parameter ID and value and stores in global parameter dictionary.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The parameter XML element

    Notes
    -----
    Adds the parameter to the global param dictionary.
    """
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
    """Process a single reaction element from SBML.

    Extracts reaction attributes including name, ID, flux bounds, and reversibility.
    Processes child elements for reactants, products, and gene associations.
    Applies filtering logic for exchange and demand reactions based on global settings.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The reaction XML element

    Raises
    ------
    SystemExit
        If a duplicate reaction ID is encountered

    Notes
    -----
    Adds the reaction to the global reactions dictionary if it passes filters.
    May skip or modify reactions based on include_ex, include_dm, and allow_dm_rev settings.
    Increments skip_ex_count, skip_dm_count, or fix_dm_count as appropriate.
    """
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
    """Process the listOfReactants element for a reaction.

    Iterates through speciesReference child elements and adds each reactant
    to the reaction with negative stoichiometric coefficients.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The listOfReactants XML element
    react : Reaction
        The reaction object to which reactants are added
    """
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "      Tag " + tag)
        if tag == 'speciesReference':
            proc_speciesReference(elt, react, -1.0)

def proc_listOfProducts (head, react):
    """Process the listOfProducts element for a reaction.

    Iterates through speciesReference child elements and adds each product
    to the reaction with positive stoichiometric coefficients.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The listOfProducts XML element
    react : Reaction
        The reaction object to which products are added
    """
    for elt in head:
        tag = rmbrace(elt.tag)
        debug (3, "      Tag " + tag)
        if tag == 'speciesReference':
            proc_speciesReference(elt, react, +1.0)

def proc_speciesReference (head, react, sign):
    """Process a single speciesReference element (reactant or product).

    Extracts species ID and stoichiometry, then adds to the reaction with the
    appropriate sign (negative for reactants, positive for products).

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The speciesReference XML element
    react : Reaction
        The reaction object to which the species is added
    sign : float
        Sign multiplier (-1.0 for reactants, +1.0 for products)

    Raises
    ------
    SystemExit
        If the referenced species ID is not found in the global species dictionary
    """
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
    """Process a geneProductAssociation element for a reaction.

    Recursively processes 'or', 'and', and 'geneProductRef' child elements.
    If any gene product reference is found, marks the reaction as gene-indicated.

    Parameters
    ----------
    head : xml.etree.ElementTree.Element
        The geneProductAssociation, or, or and XML element
    react : Reaction
        The reaction object to mark as gene-indicated

    Notes
    -----
    This implementation marks the reaction as gene-indicated if any gene product
    exists, without tracking the full logical structure of or/and relationships.
    """
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
    """Write the parsed metabolic model to a PLD format file.

    Generates a PLD file with MET (metabolite) and REACTION entries suitable for
    flux balance analysis and gap-filling tools. Reversible reactions are output
    as separate forward and reverse entries with negated coefficients.

    Parameters
    ----------
    fn : str
        Output filename for the PLD file
    xml_fn : str
        Input XML filename (for header comment)

    Notes
    -----
    File format:
    - Header comments with timestamp and conversion settings
    - MET lines: 'MET <id> "<name>"'
    - REACTION lines: 'REACTION <id> <obj_coeff> <upper_bound> <species coeffs> [__fullname__ "<name>"]'
    - Reversible reactions generate additional '<id>__rev' entries

    Optionally strips 'R_' prefix from reaction IDs and 'M_' prefix from metabolite IDs
    if fix_name is True.
    """
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
    """Print error message to stderr and exit with error code 1.

    Parameters
    ----------
    message : str
        Error message to display
    """
    print(message, file=sys.stderr)
    exit (1)

def usage (msg = ""):
    """Print usage information and exit with error code 1.

    Displays command-line syntax, argument descriptions, and switch options.

    Parameters
    ----------
    msg : str, optional
        Additional message to display before usage information (default: '')

    Notes
    -----
    Command-line switches:
    - -j #: Set default reaction objective coefficient
    - -g #: Set gene-indicated reaction objective coefficient
    - -rm: Don't remove R_ and M_ prefixes from names
    - -ex: Exclude exchange reactions
    - -dm: Exclude demand reactions
    - +dmrev: Allow reversible demand reactions
    """
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
    """Main entry point for SBML to PLD conversion.

    Parses command-line arguments, reads and processes the SBML XML file,
    applies objective coefficients, and writes the output PLD file.
    Prints summary statistics about skipped and modified reactions.

    Raises
    ------
    SystemExit
        If required arguments are missing, invalid arguments provided, or
        processing errors occur (duplicate IDs, missing objective reactions, etc.)

    Notes
    -----
    Command-line usage:
        sbmlcvt.py [-j #] [-g #] [-rm] [-ex] [-dm] [+dmrev] file.xml [out_fn]

    See usage() function for detailed argument descriptions.
    """
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
