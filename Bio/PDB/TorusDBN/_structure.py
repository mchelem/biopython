# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

""" Methods to deal with amino acids structure conversions. 

The secondary structure labels are listed in the paper:
A series of PDB related databases for everyday needs.
    Joosten RP, Te Beek TAH, Krieger E, Hekkelman ML, Hooft RWW,
    Schneider R, Sander C, Vriend G, NAR 2010; doi: 10.1093/nar/gkq1105.

"""

import inspect
import numpy

from Bio.PDB import StructureBuilder
from Bio.PDB.Polypeptide import PPBuilder, one_to_index, index_to_one, \
    index_to_three

from Bio.PDB.TorusDBN.TorusDBNExceptions import TorusDBNException



_DSSP_TO_INDEX = {'H':0, 'G':0, 'I':0, 'E':1, 'B':1}
_INDEX_TO_DSSP = {0: 'H', 1:'E', 2:'C'}


def index_to_dssp(index):
    """ Translate a secondary structure index into label.
        - H = Helix 
        - E = Strand
        - C = Coil.     
    
    >>> index_to_dssp(0)
    'H'
    
    @param index: Index of the secondary structure
    @type index: int
    
    @rtype: str
    @return: Secondary structure label.
    """
    try:
        return _INDEX_TO_DSSP.get(index)
    except ValueError:
        raise TorusDBNException(
            "Could not translate index %d to secondary structure label." %
            index
        )
    
    
def dssp_to_index(ss_label):
    """ Translate a secondary structure label into an index.
        - Helix
            - H = Alpha helix (4-12)
            - G = 3-10 helix
            - I = pi helix   
        - Strand
            - E = Strand
            - B = Isolated beta-bridge residue
        - Coil       
    
    >>> dssp_to_index('H')
    0
    
    @param ss_label: Label of the secondary structure.
    @type ss_label: str
    
    @rtype: int
    @return: Index of the secondary structure.
    """
    return _DSSP_TO_INDEX.get(ss_label, 2)


def aa_sequence_to_indexes(aa_seq):
    """ Convert a sequence of amino acids to a list of integers.
    
    @param aa_seq: A sequence of amino acids.
    @type aa_seq: str, list(str) or list(int)
    
    @rtype: list(int)
    @return: list of amino acid codes.
    """
    aa_indexes = []
    if len(aa_seq) > 0:
        if isinstance(aa_seq, list):
            # Already a list of integers
            if isinstance(aa_seq[0], int):
                aa_indexes =  aa_seq
            # A list of strings
            elif isinstance(aa_seq[0], str):
                aa_indexes = [one_to_index(aa) for aa in aa_seq]
        elif isinstance(aa_seq, str):
            aa_indexes = [one_to_index(aa) for aa in aa_seq]
        else:
            raise TorusDBNException(
                "Could not translate amino acid sequence into amino acids indexes.")
    return aa_indexes


def ss_sequence_to_indexes(ss_seq):
    """ Convert a sequence of secondary structures to a list of integers.
    
    @param ss_seq: A sequence of secondary structures.
    @type ss_seq: str, list(str) or list(int)
    
    @rtype: list(int)
    @return: list of secondary structure codes.
    """
    ss_indexes = []
    if len(ss_seq) > 0:
        if isinstance(ss_seq, list):
            if isinstance(ss_seq[0], int):
                ss_indexes = ss_seq
            elif isinstance(ss_seq[0], str):
                ss_indexes = [dssp_to_index(ss) for ss in ss_seq]
        elif isinstance(ss_seq, str):
            ss_indexes = [dssp_to_index(ss) for ss in ss_seq]
        else:
            raise TorusDBNException(
                "Could not translate amino acid secondary structure into indexes.")
    return ss_indexes
    
    
def indexes_to_aa_sequence(aa_indexes):
    """ Construct an amino acid sequence from amino acid codes.
    
    @param aa_indexes: amino acid codes.
    @type aa_indexes: list
    
    @rtype: list
    @return: amino acids identified by the one letter code.
    
    """
    output = [index_to_one(aa_index) for aa_index in aa_indexes]
    return ''.join(output)


def indexes_to_ss_sequence(ss_indexes):
    """ Construct a secondary structure sequence from structure indexes.
    
    @param ss_indexes: Secondary structure indexes.
    @type ss_indexes: list
    
    @rtype: list
    @return: secondary structure identifiers (H: Helix, E: Strand, C: Coil).
    
    """
    output = [index_to_dssp(ss_index) for ss_index in ss_indexes]
    return ''.join(output)
    
    
def build_structure(coords, aa_seq, atom_names, superimpose_reference=None):
    """ Build a protein structure object from amino acid sequence and 
        coordinates.
    
    @param coords: Coordinates for the atoms in the structure.
    @type coords: list(float)
    
    @param aa_seq: Sequence of amino acids.
    @type aa_seq: list(str)
    
    @param atom_names: The name of the atoms to be used.
    @type atom_names: list(str)
    
    @param superimpose_reference: Superimpose structure on another.
    @type superimpose_reference: Structure
    
    @rtype: Structure
    @return: The protein structure object.
    
    """
    sb = StructureBuilder.StructureBuilder()
    sb.init_structure('0')
    sb.init_model(0)
    sb.init_chain('0')
    sb.init_seg('')
    
    for i in xrange(len(coords)):
        res_index = i // len(atom_names)
        if (i % len(atom_names)) == 0:
            sb.init_residue(index_to_three(
                aa_seq[res_index]), ' ', res_index+1, ' ')
        
        coord = coords[i]
        name_index = i % len(atom_names)

        if not (index_to_three(aa_seq[res_index]) == 'GLY' and atom_names[name_index] == 'CB'):
            # Compensate for strange bug in newer versions of biopython (>= 1.53)
            if 'element' in inspect.getargspec(sb.init_atom)[0]:
                sb.init_atom(atom_names[name_index], numpy.array(coord), 0.0, 0.0, ' ', 
                    atom_names[name_index], None, atom_names[name_index][0])
            else:
                sb.init_atom(atom_names[name_index], array(coord), 0.0, 0.0, ' ',
                    atom_names[name_index])
                    
    structure = sb.get_structure()
    if superimpose_reference:
        structure = _superimpose_structures(structure, superimpose_reference)
    
    return structure
    
    
    def _superimpose_structures(structure, superimpose_reference):
        """ Superimpose two structures.
        @param structure: Structure to be superimposed.
        @type structure: Structure
        
        @param superimpose_reference: Structure to be superimposed.
        @type superimpose_reference: Structure
        
        @rtype: Structure
        @return: The superimposed structure.
        
        """        
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
        
