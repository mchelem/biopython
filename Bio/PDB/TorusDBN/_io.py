# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import os
import math
import numpy

from mocapy.framework import eMISMASK

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Residue import Residue
from Bio.PDB.Polypeptide import PPBuilder, three_to_index
from Bio.PDB.Vector import calc_dihedral

from Bio.PDB.TorusDBN.TorusDBNExceptions import TorusDBNException, \
    TorusDBNChainBreakException
from Bio.PDB.TorusDBN._structure import dssp_to_index



def create_sequence_from_file(chain_pdb, missing_residues, quiet_parser=True):
    """ Read a PDB file and creates a sequence and mismask to represent 
    its content.
    
    @param chain_pdb: The PDB file to be read.
    @type chain_pd: str
    
    @param missing_residues: A dictionary with the missing residues.
    @type missing_residues: dict
    
    @param quiet_parser: Disable PDBParser warnings.
    @type quiet_parser: bool
        
    """
    sequences = []
    mismasks = []
            
    output_data = []
    output_mismask = []
    
    parser = PDBParser(QUIET=quiet_parser)
    structure = parser.get_structure("X", chain_pdb)

    dssp = DSSP(model=structure[0], pdb_file=chain_pdb)        

    # Loop over residues in peptides
    ppb = PPBuilder()
    pp_list = ppb.build_peptides(structure[0])
    
    if len(pp_list) == 0:
        raise TorusDBNBuildPolypeptideException(
            "Could not create a list of Polypeptide objects from the file %s." 
            % (chain_pdb)
        )
    else:
       pp_chains, chain_missing_residues = _get_pp_with_chain_break(
            chain_pdb, pp_list, missing_residues)
                   
    for pp_index, pp in enumerate(pp_chains):
        phi_psi_list = pp.get_phi_psi_list()
        missing_residues = chain_missing_residues[pp_index]
                    
        for i in xrange(1, len(phi_psi_list)-1):
            seq = [0] * 6
            mism = [eMISMASK.MOCAPY_HIDDEN] + 4 * [eMISMASK.MOCAPY_OBSERVED]
            
            # Amino acid
            res = pp[i]
            res_name = res.get_resname()
            res_index = res.get_id()[1]

            aa_index = three_to_index(res_name)
            
            if res_index in missing_residues:
                seq[3] = aa_index
                mism[1] = eMISMASK.MOCAPY_MISSING # angles unknown
                mism[3] = eMISMASK.MOCAPY_MISSING # ss unknown
                mism[4] = eMISMASK.MOCAPY_MISSING # cis unknown  
            else:
                # Secondary Structure
                ss = res.xtra["SS_DSSP"]
                ss_index = dssp_to_index(ss)                
                                 
                # Angles
                if None in phi_psi_list[i]:
                # Previous or next residue missing, therefore angles are
                # Unknown
                    mism[1] = eMISMASK.MOCAPY_MISSING                                                 
                else:
                    seq[1:3] = phi_psi_list[i]
                    
                seq[3] = aa_index
                seq[4] = ss_index
                
                 # Cis/Trans information   
                if (res_index - 1) in missing_residues:
                    mism[4] = eMISMASK.MOCAPY_MISSING # cis unknown   
                else:                                        
                    seq[5] = _get_peptide_bond_conformation(res, pp[i-1])

            output_data.append(seq)
            output_mismask.append(mism)

    sequences.append(numpy.array(output_data))
    mismasks.append(numpy.array(output_mismask, dtype = numpy.uint))
    return sequences, mismasks    


def read_missing_residues(filename):  
    """ Read missing residues from a file.
    
    The file format is the following:
    # file  chain_id    position    residue
    'mychain.pdb'  1           23         ALA   # one position
    'mychain.pdb'  3          2:4         LEU   # position range
    'other.pdb     1          1,3         GLY   # position list
    
    @param filename: The file describing the missing residues.
    @type filename: str
    
    @rtype: dict
    @return: Dictionary of residues whose key is the (file, chain_id, position).
    
    """
    def get_chain_index(chain_index):
        if chain_index.find(':'):
            index_range = chain_index.split(':') 
            chain_index_list = range(
                int(index_range[0]), 
                int(index_range[-1]) + 1
            )
        else:
            chain_index_list = [int(index)]
        return chain_index_list
            
    try:
        missing_residues = {}
        if not os.path.isfile(filename):
            raise TorusDBNException(
                "Could not open file %s for reading missing residues." % filename)
       
        residues_file = open(filename)
        for i, line in enumerate(residues_file.readlines()):
            line = line.strip()
            if not line.startswith("#") and len(line) > 1:
                try: 
                    chain_pdb, chain_id, res_index, res_name = \
                        line.split()
                        
                    chain_index_list = []
                    
                    if chain_id.find(',') != -1:
                        chain_list = chain_id.split(',') 
                        for chain in chain_list:
                            chain_index_list += get_chain_index(chain)                                
                    else:
                        chain_index_list += get_chain_index(chain_id)
                    
                    for chain_index in chain_index_list:
                        missing_residues[
                            (chain_pdb, int(chain_index), int(res_index))
                        ] = Residue(id=('', int(res_index), ''), resname=res_name, segid='')
                except ValueError:
                    TorusDBNException(
                        "Could not read missing residues file %s at line %d." % 
                        (filename, i + 1)
                    )
    except IOError:
        TorusDBNException("Could not open file %s for reading missing residues." % filename)
    return missing_residues
        
        
def _get_missing_residue(missing_residues, chain_pdb, chain_index, residue_index):  
    """ Get the missing residues corresponding to a file, chain and position.
    
    @param missing_residues: Dictionary of residues indexed by 
        (file, chain_id, position).
    @type missing_residues: dict
    
    @param chain_pdb: The filename where the chain is described.
    @type chain_pdb: str
    
    @param chain_index: The index of the chain in the file.
    @type chain_index: int
    
    @param residue_index: The position of the residue in the chain
    @type residue_index: int
    
    @rtype: str
    @return: The missing residue three letter identifier.
    """
    chain_pdb = os.path.split(chain_pdb)[-1]
    try:
        return missing_residues[(chain_pdb, chain_index, residue_index)]
    except:
        raise TorusDBNChainBreakException(
            "Chain break in file %s, chain %d, residue %d, could not be handled." % 
            (chain_pdb, chain_index, residue_index)
        )
    
    
def _get_pp_with_chain_break(chain_pdb, pp_list, missing_residues_dict):
    """ Get the polypeptides chains for a given chain file, adding missing
    residues.
    
    @param chain_pdb: The file containing the chains.
    @type chain_pdb: str
    
    @param pp_list: List of the polypeptide chains.
    @type pp_list: list(Polypeptide)
    
    @param missing_residues_dict: Dictionary of residues indexed by 
        (file, chain_id, position).
    @type missing_residues_dict: dict
    
    @rtype: tuple(list(Polypeptide), list(int))
    @return: Polypeptide chains and list with missing residues position.
    """
    pp_chains = []
    missing_residues = []
                   
    for pp in pp_list:
        if pp[0].get_id()[1] == 1:
            pp_chains.append(pp)  
            missing_residues.append([])                                  
        else:
            last_residue = pp_chains[-1][-1].get_id()[1] + 1
            first_residue = pp[0].get_id()[1]                    
            missing_residues_index = range(last_residue, first_residue)
            missing_residues[-1] += missing_residues_index
            
            for res_index in missing_residues_index:
                pp_chains[-1].append(_get_missing_residue(
                    missing_residues_dict, chain_pdb, len(pp_chains), res_index))
            pp_chains[-1] += pp

    return pp_chains, missing_residues


def _get_peptide_bond_conformation(res, prev_res):
    """ Get the conformation of the peptide bond.
    
    @param res: A residue
    @type res: Residue
    
    @param res: The previous residue
    @type res: Residue
    
    @rtype: int (cis: 0, trans: 1)
    @return: Conformation of the peptide bond.
    
    """
    CIS = 0
    TRANS = 1
    
    CA_prev = prev_res['CA'].get_vector()
    C_prev = prev_res['C'].get_vector()
    N = res['N'].get_vector()
    CA = res['CA'].get_vector()
    dihedral = calc_dihedral(CA_prev, C_prev, N, CA)

    if abs(dihedral) < math.pi/4:
        return CIS
    else:
        return TRANS
