# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from time import time
import numpy

from mocapy.framework import DBN, NodeFactory, mocapy_seed
from mocapy.inference import SampleInfEngineHMM, LikelihoodInfEngineHMM

from Bio.PDB.PDBIO import PDBIO

from Bio.PDB.TorusDBN._structure import aa_sequence_to_indexes, ss_sequence_to_indexes, \
    indexes_to_aa_sequence, indexes_to_ss_sequence, build_structure
from Bio.PDB.TorusDBN._geometry import calculate_coordinates


# Conformation of the peptide bond (fixed at 180 (trans) or 0 (cis))
CIS = 0
TRANS = 1


class TorusDBNModel(object):
    """ TorusDBN is a probabilistic model that attempts to efficiently explore the
    conformational space in a manner that reflects its stability. The angular 
    degrees of freedom are restricted to values found in backbones of native 
    structures.
    
    """

    def __init__(self, seed=int(time())):    
        """
        @param seed: Seed for the random number generator.
        @type seed: int
        
        """
        # Mocapy config
        mocapy_seed(seed)   
              
        # Node sizes
        self.size_h = 55 # Default hidden node size     
        self.size_aa = 20 # The 20 amino acids       
        self.size_ss = 3 # helix, strand, coil                
        self.size_cis = 2 # cis, trans

        # Default Model
        self.dbn = self._create_dbn()
        
        # Sampling sequence and mismask
        self.sequence = None
        self.mismask_sample = None
        self.mismask = None
        self._must_update_sequence = False
        
        # Model parameters
        self._aa = []
        self._ss = []
        self._cis = []
        self._angles = []            
                
     
    def load_dbn(self, filename):
        """ Load a dynamic Bayesian network from file.
        
        @param filename: The bayesian network file.
        @type filename: str    
            
        """
        self.dbn = DBN()
        self.dbn.load(filename)
        
        
    def save_dbn(self, filename):
        """ Save a dynamic Bayesian network to a file.
        
        @param filename: File to which the Bayesian network will be written.
        @type filename: str
                
        """
        self.dbn.save(filename) 
        
    
    def create_dbn(self, hidden_node_size):
        """ Create a new dynamic Bayesian network.
        
        @param hidden_node_size: Size of the hidden node.
        @type hidden_node_size: int        
        """
        self.size_h = hidden_node_size
        self.dbn = self._create_dbn()        
        
        
    def _create_dbn(self):
        """ Create the TorusDBN dynamic Bayesian network.
                   .-------.                         .-------.
                   |  h1   |------------------------>|  h2   |------------>(...)
           +-------'-------'-------+         +-------'-------'-------+
           |       |       |       |         |       |       |       |
           |       |       |       |         |       |       |       |
        .--v--. .--v--. .--v--. .--v--.   .--v--. .--v--. .--v--. .--v--.
        |  d  | |  a  | |  s  | |  c  |   |  d  | |  a  | |  s  | |  c  |
        '-----' '-----' '-----' '-----'   '-----' '-----' '-----' '-----'
                    Slice 1                            Slice 2
       
        h: hidden node
        d: dihedral angles (-180 to 180, a point on the torus)
        a: amino acids
        s: secondary structure (helix, strand, coil)
        c: conformation of the peptide bond (assumed to be fixed at 180 
            (trans) or 0 (cis))
            
        """
        # Nodes in slice 1
        hidden_1 = NodeFactory.new_discrete_node(self.size_h, "h1")        
        mus = numpy.zeros((self.size_h, 2))
        kappas = numpy.ones((self.size_h, 3)) * 0.5        
        dihedral_angles = NodeFactory.new_vonmises2d_node("d", mus, kappas)        
        amino_acids = NodeFactory.new_discrete_node(self.size_aa, "a")
        secondary_structure = NodeFactory.new_discrete_node(self.size_ss, "s")
        cis = NodeFactory.new_discrete_node(self.size_cis, "c")
     
        # Nodes in slice 2
        hidden_2 = NodeFactory.new_discrete_node(self.size_h, "h2")

        dbn = DBN()        
        dbn.set_slices(
            [hidden_1, dihedral_angles, amino_acids, secondary_structure, cis],
            [hidden_2, dihedral_angles, amino_acids, secondary_structure, cis],
        )        
        dbn.add_inter("h1", "h2")
        dbn.add_intra("h1", "d")
        dbn.add_intra("h1", "a")
        dbn.add_intra("h1", "s")
        dbn.add_intra("h1", "c")

        dbn.construct()
        return dbn
        
        
    def get_log_likelihood(self):
        """
        Calculate the log likelihood.
        
        @rtype: float
        @return: Log likelihood
        """        
        if self.sequence is not None:
            hmm_ll_calculator = LikelihoodInfEngineHMM(
                dbn=self.dbn, hidden_node_index=0, check_dbn=False)
            ll = hmm_ll_calculator.calc_ll(self.sequence, self.mismask)
        else:
            raise TorusDBNException(
                "TorusDBN.get_log_likelihood() can only be "
                "called after a sample has been drawn.")            
        return ll
        
        
    def sample(self):
        """ Generate a new sample. """
        if self._must_update_sequence:
            self._update_sequence()            
        
        inf_engine = SampleInfEngineHMM(
            self.dbn, self.sequence, self.mismask_sample, hidden_node_index=0)                        
        self.sample_data = inf_engine.sample_next()     
        
        self._save_sample(self.sample_data)   
                
        
    def _save_sample(self, sample):
        """ Save the generated sample as model parameters. 
        
        @param sample: Sample to be saved
        @type sample: numpy.array
        
        """ 
        self._angles = numpy.array([sample[:,1], sample[:,2]]).transpose()
        self._aa = sample[:,3]
        self._ss = sample[:,4]
        self._cis = sample[:,5]       

        
    def _update_sequence(self):  
        """ Update the sequence with the current model parameters. """
        sequence_len = max(len(self._aa), len(self._ss), len(self._cis), 
            len(self._angles))   
        num_nodes = 5

        self.sequence = numpy.zeros((sequence_len, num_nodes + 1)) 
        self.mismask_sample = numpy.ones((sequence_len, num_nodes), dtype=numpy.uint)
        self.mismask = numpy.zeros((sequence_len, num_nodes), dtype=numpy.uint)
        self.mismask[:,0] = 1
                
        if len(self._angles) > 0:    
            angles = self._angles.transpose()
            self.mismask_sample[:,1] = numpy.zeros(sequence_len) 
            self.sequence[:,1] = numpy.array(angles[0])  
            self.sequence[:,2] = numpy.array(angles[1])                       

        self._set_sequence_and_mismask(self._aa, 2, sequence_len) 
        self._set_sequence_and_mismask(self._ss, 3, sequence_len) 
        self._set_sequence_and_mismask(self._cis, 4, sequence_len) 
              
        self._must_update_sequence = False 
        
        
    def _set_sequence_and_mismask(self, parameter, position, sequence_len):
        """ Create a sequence and mismask for a given model parameter. 
        
        @param parameter: Model parameter
        @type parameter: list(Number)
        
        @param position: Position of the parameter in the mismask array.
        @type position: int
        
        @param sequence_len: Size of the sequence to be created.
        @type sequence_len: int
        
        """
        if len(parameter) > 0:
            self.mismask_sample[:,position] = numpy.zeros(
                sequence_len, dtype=numpy.uint) 
            self.sequence[:,position + 1] = numpy.array(parameter)        

    
    def save_structure(self, filename, superimpose_structure=None):  
        """ Save the generated sample to a pdb file.
        
        @param filename: The file to which the structure will be written.
        @type filename: str
        
        @param superimpose_structure: Structure to which the sampled structure
            is to be superimposed.
        @type superimpose_structure: Structure
        
        """
        io = PDBIO()
        io.set_structure(self.get_structure())
        io.save(filename)
        
    
    # -- Getters and setters for the model parameters
    
    def set_aa(self, aa):
        """
        @param aa: Amino acid sequence.
        @type aa: str, list(str), list(int)
        
        >> model = TorusDBNModel()
        >> model.set_aa('ACDE') # amino acid one letter identifiers
        >> model.set_aa(['A', 'C', 'D', 'E']
        >> model.set_aa([1, 2, 3, 4]) # amino acid codes
        
        """
        self._aa = aa_sequence_to_indexes(aa)
        self._must_update_sequence = True
        
        
    def set_ss(self, ss):
        """
        @param ss: Secondary structure of the sequence.
        @type ss: str, list(str), list(int)
        
        >> model = TorusDBNModel()
        >> model.set_ss('ECCH') # DSSP identifiers
        >> model.set_aa(['E', 'C', 'C', 'H']
        >> model.set_aa([1, 2, 2, 0]) # DSSP codes
        
        """
        self._ss = ss_sequence_to_indexes(ss)
        self._must_update_sequence = True
        
        
    def set_cis(self, cis):
        """
        @param cis: Conformation of the peptide bond (trans or cis)

        @type cis: int (trans: 0, cis: 1)
        
        >> model = TorusDBNModel()
        >> model.set_cis([1, 1, 0, 0]) 
                
        """
        
        self._cis = cis
        self._must_update_sequence = True
        
    
    def set_angles(self, angles):
        """
        @param angles: Dihedral angles.
        @type param: list(list(float))
        
        >> model = TorusDBNModel()
        >>  model.set_aa('ACDE')
        >> model.set_angles([
            [1.49, -1.64],
            [-0.64, 1.34],
            [1.32, -0.79],
            [-0.66, 1.29],
        ])    
        
        """
        self._angles = numpy.array(angles)
        self._must_update_sequence = True
        
        
    def get_ss(self):
        """
        @rtype: str
        @return: Secondary structure labels.
        
        >> model = TorusDBNModel()
        >> model.set_aa('ACDEFGHIK')                     
        >> model.sample()
        >> model.get_ss()
        'ECCE'    
            
        """
        return indexes_to_ss_sequence(self._ss)
        
    
    def get_aa(self):
        """
        @rtype: str
        @return: Amino acids represented in one letter format.
        
        >> model = TorusDBNModel()
        >> model.set_aa('ACDEFGHIK')                     
        >> model.sample()
        >> model.get_aa()
        'ACDEFGHIK'    
            
        """
        return indexes_to_aa_sequence(self._aa)
        
    
    def get_angles(self):
        """
        @rtype: numpy.array
        @return: Dihedral angles
        
        >> model = TorusDBNModel()
        >> model.set_aa('ACDE')                     
        >> model.sample()
        >> model.get_angles()
        [[ 1.49 -1.64]
         [-0.64  1.34]
         [ 1.32 -0.79]
         [-0.66  1.29]]
         
        """
        return self._angles
        
    
    def get_cis(self):
        """
        @rtype: numpy.array
        @return: Peptide bond conformation (trans: 0, cis: 1)
        
        >> model = TorusDBNModel()
        >> model.set_aa('ACDE')                     
        >> model.sample()
        >> model.get_cis()
        [ 1.  0.  1.  0.]
        
        """
        return self._cis
    
        
    def get_structure(self, superimpose_structure=None):
        """
        @rtype: Structure
        @return: Structure that represents the aminoacid sequence
        """
        return build_structure(
            self.get_coordinates(), 
            self._aa, 
            ['N', 'CA', 'C', 'CB', 'O'],                                  
            superimpose_structure,
        )
        
        
    def get_coordinates(self):
        """ Get sample coordinates according to the amino acid residues, 
        angles and cis/trans conformation.
    
        @rtype: numpy.array
        @return: The calculated coordinates.
        """
        return calculate_coordinates(self._angles, self._cis, self._aa)
        
