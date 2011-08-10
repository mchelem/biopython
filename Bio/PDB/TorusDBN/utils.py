# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from Bio.PDB.Polypeptide import PPBuilder, one_to_index, index_to_one

__dssp_to_index = {'H':0, 'G':0, 'I':0, 'E':1, 'B':1}
__index_to_dssp = {0: 'H', 1:'E', 2:'C'}

def index_to_dssp(index):
    """ 
    Translate index to secondary structure label
    H = Helix, 
    E = Strand,
    C = Coil.     
    
    >>> index_to_dssp(0)
    'H'
    """
    return __index_to_dssp.get(index)
    
    
def dssp_to_index(ss):
    """ 
    Translate secondary structure label to indes
    H = Alpha helix (4-12), 
    G = 3-10 helix,
    I= = pi helix,
    
    E = Strand,
    B = Isolated beta-bridge residue.
    
    >>> dssp_to_index('H')
    0
    """
    return __dssp_to_index.get(ss, 2)


def aa_to_list(aa_seq):
    '''Convert aa_seq to list of integers if not already so'''
    if len(aa_seq) == 0:
        return aa_seq
    new_aa_seq = aa_seq
    if isinstance(aa_seq, list):
        # Determine type based on first element
        if isinstance(aa_seq[0], int):
            pass
        elif isinstance(aa_seq[0], str):
            new_aa_seq = [one_to_index(aa) for aa in aa_seq]
    elif isinstance(aa_seq, str):
        new_aa_seq = [one_to_index(aa) for aa in aa_seq]
    return new_aa_seq


def ss_to_list(ss_seq):
    '''Convert ss_seq to list of integers if not already so'''
    if len(ss_seq) == 0:
        return ss_seq
    new_ss_seq = ss_seq
    if isinstance(ss_seq, list):
        # Determine type based on first element
        if isinstance(ss_seq[0], int):
            pass
        elif isinstance(ss_seq[0], str):
            new_ss_seq = [dssp_to_index(ss) for ss in ss_seq]
    elif isinstance(ss_seq, str):
        new_ss_seq = [dssp_to_index(ss) for ss in ss_seq]
    return new_ss_seq
