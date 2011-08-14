# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from Bio.PDB.Polypeptide import PPBuilder, one_to_index, index_to_one
from Bio.PDB import StructureBuilder
from Bio.PDB.Polypeptide import index_to_three
import inspect
import numpy

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
    
    
def build_structure(coords, aa_seq, atom_names, superimpose_reference=None):
    """
    Build protein structure object    
    """
    sb = StructureBuilder.StructureBuilder()
    sb.init_structure('0')
    sb.init_model(0)
    sb.init_chain('0')
    sb.init_seg('')
    for i in xrange(len(coords)):
        res_index = i // len(atom_names)
        if i%len(atom_names)==0:
                sb.init_residue(index_to_three(aa_seq[res_index]), ' ', res_index+1, ' ')
        coord = coords[i]
        name_index = i%len(atom_names)

        if not (index_to_three(aa_seq[res_index]) == 'GLY' and atom_names[name_index] == 'CB'):
            # Compensate for strange bug in newer versions of biopython (>=1.53)
            if 'element' in inspect.getargspec(sb.init_atom)[0]:
                sb.init_atom(atom_names[name_index], numpy.array(coord), 0.0, 0.0, ' ', 
                    atom_names[name_index], None, atom_names[name_index][0])
            else:
                sb.init_atom(atom_names[name_index], array(coord), 0.0, 0.0, ' ',
                    atom_names[name_index])
    structure = sb.get_structure()

    if superimpose_reference:
        ppb = PPBuilder()
        sup = Superimposer()
        pp_reference = ppb.build_peptides(superimpose_reference)[0]
        pp_structure = ppb.build_peptides(structure)[0]

        # CA only
        fixed = pp_reference.get_ca_list()
        moving = pp_structure.get_ca_list()
        moving_all = Selection.unfold_entities(structure, "A")        
        sup.set_atoms(fixed, moving)
        sup.apply(moving_all)
    
    return structure
    
    
def build_sequence_aa(aa_seq):
    """
    Construct an amino acid sequence
    """
    output = [index_to_one(aa) for aa in aa_seq]
    return ''.join(output)


def build_sequence_ss(ss_seq):
    """
    Construct a secondary structure sequence
    """
    output = [index_to_dssp(ss) for ss in ss_seq]
    return ''.join(output)
