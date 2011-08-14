# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from mocapy.framework import DBN, NodeFactory, mocapy_seed
from mocapy.inference import SampleInfEngineHMM, LikelihoodInfEngineHMM

from Bio.PDB.PDBIO import PDBIO

from Bio.PDB.TorusDBN._utils import aa_to_list, ss_to_list
from Bio.PDB.TorusDBN._utils import build_structure, build_sequence_aa, \
    build_sequence_ss
from Bio.PDB.TorusDBN._geometry import get_coordinates_from_angles


import time
import numpy



class TorusDBNModel(object):
    """
    TorusDBN is a probabilistic model that attempts to efficiently explore the
    conformational space in a manner that reflects its stability. The angular 
    degrees of freedom are restricted to values found in backbones of native 
    structures.
    """

    def __init__(self, seed=int(time.time())):    
        # Mocapy config
        mocapy_seed(seed)   
              
        # Node sizes
        self.size_h = 55 # Default hidden node size     
        self.size_aa = 20 # The 20 amino acids       
        self.size_ss = 3 # helix, strand, coil                
        self.size_cis = 2 # cis, trans

        # Default Model
        self.dbn = self.__create_dbn()
        
        # Sampling sequence and mismask
        self.sequence = None
        self.mismask_sample = None
        self.mismask = None
        self.__must_update_sequence = False
        
        # Model parameters
        self.__aa = []
        self.__ss = []
        self.__cis = []
        self.__angles = []            
                
     
    def load_dbn(self, filename):
        self.dbn = DBN()
        self.dbn.load(filename)
        
        
    def save_dbn(self, filename):
        self.dbn.save(filename) 
        
    
    def create_dbn(self, hidden_node_size):
        self.size_h = hidden_node_size
        self.dbn = self.__create_dbn()        
        
        
    def __create_dbn(self):
        """
        Create the TorusDBN dynamic Bayesian network:
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
        if self.sequence is not None:
            hmm_ll_calculator = LikelihoodInfEngineHMM(
                dbn=self.dbn, hidden_node_index=0, check_dbn=False)
            ll = hmm_ll_calculator.calc_ll(self.sequence, self.mismask)
        else:
            raise TorusDBNException("TorusDBN.get_log_likelihood() can only be called after a sample has been drawn.")
            
        return ll
        
        
    def sample(self):
        
        if self.__must_update_sequence:
            self.__update_sequence()            
        
        inf_engine = SampleInfEngineHMM(
            self.dbn, self.sequence, self.mismask_sample, hidden_node_index=0)                        
        self.sample_data = inf_engine.sample_next()     
        
        self.__save_sample(self.sample_data)   
        
        return self.sample_data
        
        
    def __save_sample(self, sample):
        self.__angles = numpy.array([sample[:,1], sample[:,2]]).transpose()
        self.__aa = sample[:,3]
        self.__ss = sample[:,4]
        self.__cis = sample[:,5]       

        
    def __update_sequence(self):  
        """
        Sequence:  hidden phi psi aa ss cis (0 - 5)
        Mismask:  hidden vonmises aa ss cis (0 - 4)
        
        """
        sequence_len = max(len(self.__aa), len(self.__ss), len(self.__cis), 
            len(self.__angles))   
        num_nodes = 5

        self.sequence = numpy.zeros((sequence_len, num_nodes + 1)) 
        self.mismask_sample = numpy.ones((sequence_len, num_nodes), dtype=numpy.uint)
        self.mismask = numpy.zeros((sequence_len, num_nodes), dtype=numpy.uint)
        self.mismask[:,0] = 1
        
        if len(self.__angles) > 0:            
            self.mismask_sample[:,1] = numpy.zeros(sequence_len) 
            self.sequence[:,1] = numpy.array(self.__angles[0])  
            self.sequence[:,2] = numpy.array(self.__angles[1])                       

        self.__set_sequence_and_mismask(self.__aa, 2, sequence_len) 
        self.__set_sequence_and_mismask(self.__ss, 3, sequence_len) 
        self.__set_sequence_and_mismask(self.__cis, 4, sequence_len) 
              
        self.__must_update_sequence = False 
        
        
    def __set_sequence_and_mismask(self, parameter, position, sequence_len):
        if len(parameter) > 0:
            self.mismask_sample[:,position] = numpy.zeros(
                sequence_len, dtype=numpy.uint) 
            self.sequence[:,position + 1] = numpy.array(parameter)        
    
    
    def set_aa(self, aa):
        self.__aa = aa_to_list(aa)
        self.__must_update_sequence = True
        
        
    def set_ss(self, ss):
        self.__ss = ss_to_list(ss)
        self.__must_update_sequence = True
        
        
    def set_cis(self, cis):
        self.__cis = cis
        self.__must_update_sequence = True
        
    
    def set_angles(self, angles):
        self.__angles = angles
        self.__must_update_sequence = True
        
        
    def get_ss(self):
        return build_sequence_ss(self.__ss)
        
    
    def get_aa(self):
        return build_sequence_aa(self.__aa)
        
    
    def get_angles(self):
        return self.__angles
        
    
    def get_cis(self):
        return self.__cis
    
        
    def get_structure(self, superimpose_structure=None):
        return build_structure(
            self.get_coordinates(), 
            self.__aa, 
            ['N', 'CA', 'C', 'CB', 'O'],                                  
            superimpose_structure,
        )
        
        
    def get_coordinates(self):
        return get_coordinates_from_angles(self.__angles, self.__cis, self.__aa)
    
    
    def save_structure(self, filename, superimpose_structure=None):        
        io = PDBIO()
        io.set_structure(self.get_structure())
        io.save(filename)
        
