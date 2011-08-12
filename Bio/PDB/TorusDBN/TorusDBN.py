# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from mocapy.framework import DBN, NodeFactory, EMEngine, mocapy_seed
from mocapy.inference import GibbsRandom, LikelihoodInfEngineHMM, SampleInfEngineHMM

from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import PPBuilder, three_to_index, one_to_index
from Bio.PDB.Vector import calc_dihedral

from Bio.PDB.TorusDBN.utils import aa_to_list, ss_to_list, dssp_to_index
from Bio.PDB.TorusDBN.utils import build_structure, build_sequence_aa, build_sequence_ss
from Bio.PDB.TorusDBN.geometry import get_coordinates_from_angles
from Bio.PDB.TorusDBN.TorusDBNExceptions import TorusDBNBuildPolypeptideException

import numpy
import math
import time

CIS = 0
TRANS = 1



class TorusDBN(object):
    """
    TorusDBN is a probabilistic model that attempts to efficiently explore the
    conformational space in a manner that reflects its stability. The angular 
    degrees of freedom are restricted to values found in backbones of native 
    structures.
    """

    def __init__(self, seed=int(time.time())):             
        # Node sizes
        self.size_h = 55 # Default hidden node size     
        self.size_aa = 20 # The 20 amino acids       
        self.size_ss = 3 # helix, strand, coil                
        self.size_cis = 2 # cis, trans
        
        # Inference engine parameters
        self.em_steps = 100
        self.burnin_steps = 20
        self.mcmc_steps = 1
        
        # Convergence check
        self.check_convergence = False
        self.convergence_window = 80
        self.convergence_threshold = 55

        # Default Model
        self.dbn = self.__create_dbn()
        
        # Training sequences and mismasks 
        self.seq_list = []
        self.mismask_list = []
        
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
        
        # Mocapy config
        mocapy_seed(seed)
                
     
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

    
    def train(self, training_set, use_aic=False):
        self.seq_list, self.mismask_list = self.__create_sequence_and_mismask(
            training_set)      
        return self.__train(use_aic)
        
        
    def __train(self, use_aic):                         
        self.hmm_ll_calculator = LikelihoodInfEngineHMM(
            dbn=self.dbn, hidden_node_index=0, check_dbn=False)
                    
        mcmc = GibbsRandom(self.dbn)
        em = EMEngine(self.dbn, mcmc, self.seq_list, self.mismask_list, [])
        
        # EM loop
        ll_list = []
        for i in xrange(self.em_steps):
            if i == 0:
		       em.do_E_step(1, self.burnin_steps)
            else:
		       em.do_E_step(self.mcmc_steps, 0)
            
            ll_list.append(em.get_loglik())
            em.do_M_step()
            
            if self.check_convergence:
                if i > window_size:
                     if self.__get_hairiness(ll_list):
                        break
                        
        if (use_aic):
            return self.__calculate_AIC()
        else:
            return self.__calculate_BIC()

    
    def __calculate_AIC(self):
        """
        The Akaike information criterion (AIC) is a measure of the relative 
        goodness of fit of a statistical model.
        
        Akaike, Hirotugu (1974). "A new look at the statistical model 
        identification". IEEE Transactions on Automatic Control 19 (6): 
        716-723. doi:10.1109/TAC.1974.1100705. MR0423716.
        """            
        hmm_ll_calculator = LikelihoodInfEngineHMM(
            dbn=self.dbn, hidden_node_index=0, check_dbn=False)
        ll_full = hmm_ll_calculator.calc_ll(self.seq_list, self.mismask_list)        
        return 2 * ll_full - 2 * self.__get_parameter_count()
        
                
    def __calculate_BIC(self):
        """
        The Bayesian information criterion (BIC) is a criterion for model 
        selection among a finite set of models.
        
        Schwarz, Gideon E. (1978). "Estimating the dimension of a model". Annals
        of Statistics 6 (2): 461-464. doi:10.1214/aos/1176344136. MR468014.
        """            
        hmm_ll_calculator = LikelihoodInfEngineHMM(
            dbn=self.dbn, hidden_node_index=0, check_dbn=False)
        ll_full = hmm_ll_calculator.calc_ll(self.seq_list, self.mismask_list)              
        return 2 * ll_full - self.__get_parameter_count() * math.log(
            self.__get_observation_count())

    
    def __get_hairiness(self, ll_list):
        for j in xrange(self.convergence_window, len(ll_list)):
	
            match=0
            for i in xrange(self.convergence_window):               
                # Test if likelihood is oscillating
                if (ll_list[j-i] > ll_list[j-i+1] and ll_list[j-i] > ll_list[j-i-1]) \
                    or (ll_list[j-i] < ll_list[j-i-1] and ll_list[j-i] < ll_list[j-i+1]) :                
                    match += 1                
            
            if match >= self.convergence_threshold:
                return True
	
        return False	
        
    
    def __get_observation_count(self):
        observation_count = 0
        for sequence in self.seq_list:
            observation_count += sequence.shape[0]     
            
        return observation_count
        
                
    def __get_parameter_count(self):
        parameters_d = 5;
        return (self.size_h - 1) + self.size_h * (
            (self.size_h - 1) + parameters_d + (self.size_aa - 1) + 
            (self.size_ss - 1) + (self.size_cis - 1)
        )     
        
                
    def __create_sequence_and_mismask(self, training_set):
        seq_list = []
        mismask_list = []
        for filename in training_set:
            try:
                sequence, mismask = self.__get_data_from_file(filename)
                seq_list.append(sequence)
                mismask_list.append(mismask)
            except TorusDBNBuildPolypeptideException:
                print "\nCould not create polypeptide list from file %s ." % (filename)
                print "This file was not included in the training set.\n"
            
        return seq_list, mismask_list
        
    
    def __get_data_from_file(self, chain_pdb):
        output_data = []
        output_mismask = []
        
        parser = PDBParser()
        structure = parser.get_structure("X", chain_pdb)

        dssp = DSSP(model=structure[0], pdb_file=chain_pdb)

        # Loop over residues in peptides
        ppb = PPBuilder()
        pp_list = ppb.build_peptides(structure[0])

        try:
            pp = pp_list[0]
        except IndexError:
            raise TorusDBNBuildPolypeptideException(
                "Could not create a list of Polypeptide objects from the file %s." 
                % (chain_pdb)
            )

        phi_psi_list = pp.get_phi_psi_list()
        
        for i in xrange(1, len(phi_psi_list)-1):
            seq = [0] * 6

            # Amino acid
            res = pp[i]
            res_name = res.get_resname()
            aa_index = three_to_index(res_name)

            # Secondary Structure
            ss = res.xtra["SS_DSSP"]
            ss_index = dssp_to_index(ss)

            # Angles
            seq[1:3] = phi_psi_list[i]
            seq[3] = aa_index
            seq[4] = ss_index

            # Cis/Trans information:
            prev_res = pp[i-1]
            dihedral = None

            CA_prev = prev_res['CA'].get_vector()
            C_prev = prev_res['C'].get_vector()
            N = res['N'].get_vector()
            CA = res['CA'].get_vector()
            dihedral = calc_dihedral(CA_prev, C_prev, N, CA)

            if abs(dihedral) < math.pi/4:
                seq[5] = CIS
            else:
                seq[5] = TRANS

            output_data.append(seq)
            output_mismask.append([1] + 4 * [0])

        sequence = numpy.array(output_data)
        mismask = numpy.array(output_mismask, dtype = numpy.uint)
        return sequence, mismask        
            
            
    def find_optimal_model(self, training_set, use_aic=False, min_node=10, 
        max_node=90, start_int=20, end_int=5, node_samples=4, full_ll_dec=False):
            
        self.seq_list, self.mismask_list = self.__create_sequence_and_mismask(
            training_set)             
            
        max_position = 0
        start_res = start_int
        avg_full_LL = []
        
        IC_array = [[]*n for n in xrange(node_samples + 2)]
        
        # Decrease size resolution until threshold (end_int)
        while start_int >= end_int:
            # Loop over node sizes
            for i in xrange(min_node, max_node + 1, start_int):
                
                # Continues if at the maximum node size from the previous resolution 
                if (len(IC_array[0]) > 0 and i == IC_array[0][max_position]) or i <= 0:
                    continue

                # Add node-size value to header
                IC_array[0].append(i)
                IC_cum = 0
                
                if start_res == start_int:
                    avg_full_LL.append(0)
                    
                for j in xrange(1, node_samples + 1):
                    print "Training with node size = %d (sample %d)" % (i, j)
                    self.create_dbn(hidden_node_size=i)
                    IC = self.__train(use_aic)
                    IC_array[j].append(IC)
                    IC_cum += IC
                    
                    if (full_ll_dec):
                        # Save forward likelihoods in order to infer if it is decreasing
                        hmm_ll_calculator = LikelihoodInfEngineHMM(
                            dbn=self.dbn, hidden_node_index=0, check_dbn=False)
                        ll_full = hmm_ll_calculator.calc_ll(self.seq_list, self.mismask_list)
                        avg_full_LL[-1] = avg_full_LL[-1] + ll_full/self.__get_observation_count()
                    
                # Calculate mean IC for each node-size and add to array
                IC_array[node_samples + 1].append(IC_cum / node_samples)
                
                # Check if log-likelihood is decreasing 
                if (len(avg_full_LL) > 1) and (avg_full_LL[-1] < avg_full_LL[-2]) and \
                    (start_res == start_int) and full_ll_dec:
                    print "Log-likelihood is decreasing. There is no reason to test higher node sizes."
                    break
                                     
            # Column number for maximum IC value
            max_position = IC_array[node_samples + 1].index(max(IC_array[node_samples + 1])) 
            print "Optimal node size:", IC_array[0][max_position], "\n"
          
            # Update resolution
            start_int = start_int / 2
            
            # Update node limits
            min_node = IC_array[0][max_position] - start_int
            max_node = IC_array[0][max_position] + start_int
     
        IC_max_node = IC_array[0][max_position]
        
        # Final train to the optimal model
        dbn_list = []
        IC_list = []
        
        for j in xrange(node_samples):
            self.create_dbn(hidden_node_size=IC_max_node)
            IC = self.__train(use_aic)
            IC_list.append(IC)
            dbn_list.append(self.dbn)
            
        IC_max = max(IC_list)
        self.dbn = dbn_list[IC_list.index(IC_max)]
        
        print "Optimal Model:"
        print "Hidden node size =", IC_max_node
        print "IC =", IC_max
        
        return IC_max_node, IC_max
        
        
    def get_log_likelihood(self):        
        if self.sequence is not None:
            hmm_ll_calculator = LikelihoodInfEngineHMM(
                dbn=self.dbn, hidden_node_index=0, check_dbn=False)
            ll = hmm_ll_calculator.calc_ll(self.sequence, self.mismask)
        else:
            raise Exception("TorusDBN.get_log_likelihood() can only be called after a sample has been drawn.")
            
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
        
