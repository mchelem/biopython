from mocapy.framework import DBN, NodeFactory, EMEngine, mocapy_seed
from mocapy.inference import GibbsRandom, LikelihoodInfEngineHMM

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import PPBuilder, three_to_index
from Bio.PDB.Vector import calc_dihedral

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

    def __init__(self, number_of_states=55):
        """
        @param number_of_nodes: int
            The number of states for the hidden node.
        """
        # Node sizes
        self.size_h = number_of_states      
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
        
        # Model
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
        # Mocapy seed
        t = int(time.time())
        mocapy_seed(t)

        # Nodes in slice 1
        hidden_1 = NodeFactory.new_discrete_node(self.size_h, "h1")
        dihedral_angles = NodeFactory.new_vonmises2d_node("d")
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
        
        
    def train(self, training_set):
        """
        @param training_seq: list(str)
            The training set for the network.
        """        
        mcmc = GibbsRandom(self.dbn)
        self.seq_list, self.mismask_list = self.__create_sequence_and_mismask(
            training_set)
        em = EMEngine(self.dbn, mcmc, self.seq_list, self.mismask_list, [])
        #print self.dbn
        
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

    
    def calculate_AIC(self):
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
        
                
    def calculate_BIC(self):
        """
        The Bayesian information criterion (BIC) is a criterion for model 
        selection among a finite set of models.
        
        Schwarz, Gideon E. (1978). "Estimating the dimension of a model". Annals
        of Statistics 6 (2): 461-464. doi:10.1214/aos/1176344136. MR468014.
        """
        hmm_ll_calculator = LikelihoodInfEngineHMM(
            dbn=self.dbn, hidden_node_index=0, check_dbn=False)
            
        ll_full = hmm_ll_calculator.calc_ll(self.seq_list, self.mismask_list)        
        
        return 2 * ll_full - 2 * self.__get_parameter_count() * math.log(
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
            sequence, mismask = self.__get_data_from_file(filename)
            seq_list.append(sequence)
            mismask_list.append(mismask)
            
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

        pp = pp_list[0]
        phi_psi_list = pp.get_phi_psi_list()

        for i in xrange(1, len(phi_psi_list)-1):
            seq = [0] * 6

            # Amino acid
            res = pp[i]
            res_name = res.get_resname()
            aa_index = three_to_index(res_name)

            # Secondary Structure
            ss = res.xtra["SS_DSSP"]
            ss_index = self.__dssp2i(ss)

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
    
	
    def __dssp2i(self, ss):
	    "Convert SS label to index"
	    if ss in ("H", "G", "I"):
	        return 0
	    elif ss in ("E","B"):
	        return 1
	    else:
	        return 2	
