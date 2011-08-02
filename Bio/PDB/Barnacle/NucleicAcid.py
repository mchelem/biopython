# Copyright (C) 2008 by Jes Frellsen, Ida Moltke and Martin Thiim
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""
Nucleic Acid related variables, methods and classes
"""

from types import StringType

to_one_letter_code = {'  A':'A', '  C':'C', '  G':'G', '  T':'T','  U':'U'}
to_three_letter_code = {'A':'  A', 'C':'  C', 'G':'  G', 'T':'  T','U':'  U'}
purine_one_letter_codes = ['A','G']
pyrimidine_one_letter_codes = ['C','T','U']
rna_one_letter_codes = ['A','C','G','U']
rna_three_letter_codes = ['  A','  C','  G','  U']

def get_one_letter_code(residue):
    """
    Returns one letter code for a residue or converts a three letter
    code to a one letter code.
    
    @param residue: a L{Residue} object OR the three letter nucleic acid code
    @type residue: L{Residue} or string
    @return: one letter code (in UPPER)
    @rtype: L{string}

    """
    if not type(residue)==StringType:
        residue=residue.get_resname()
    residue=residue.upper()
    if to_one_letter_code.has_key(residue):
	return to_one_letter_code[residue]
    else:
	return None

def get_three_letter_code(residue):
    """
    Returns three letter code for a residue or converts a one letter
    code to a three letter code.
    
    @param residue: a L{Residue} object OR the three letter nucleic acid code
    @type residue: L{Residue} or string
    @return: three letter code (in UPPER)
    @rtype: L{string}

    """
    if not type(residue)==StringType:
        return residue.get_resname().upper()
    else:
        residue.upper()
        if to_three_letter_code.has_key(residue):
            return to_three_letter_code[residue]
        else:
            return "  X"



def is_rna_nuc(residue):
    """
    Return 1 if residue object/string is an RNA nucleotide.

    @param residue: a L{Residue} object OR the one letter nucleic acid code
    @type residue: L{Residue} or string
    @return: True/False
    @rtype: L{bool}
    """
    return (get_one_letter_code(residue) in rna_one_letter_codes)  


def is_purine(residue):
    """
    Return 1 if residue object/string is a purine

    @param residue: a L{Residue} object OR the one letter nucleic acid code
    @type residue: L{Residue} or string
    @return: True/False
    @rtype: L{bool}
    """
    return (get_one_letter_code(residue) in purine_one_letter_codes)


def is_pyrimidine(residue):
    """
    Return 1 if residue object/string is a pyrimidine

    @param residue: a L{Residue} object OR the one letter nucleic acid code
    @type residue: L{Residue} or string
    @return: True/False
    @rtype: L{bool}
    """
    return (get_one_letter_code(residue) in pyrimidine_one_letter_codes)


def has_base_type(type,residue):
    """
    Return 1 if residue object has base type 'type' 

    @param type: the one letter nucleic acid code (in UPPER)
    @type type: string
    @param residue: a L{Residue} object
    @type residue: L{Residue}
    @return: True/False
    @rtype: L{bool}
    """
    return (get_one_letter_code(residue) == type)


