# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from mocapy.framework import EMEngine, mocapy_seed, eMISMASK
from mocapy.inference import GibbsRandom, LikelihoodInfEngineHMM

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import PPBuilder, three_to_index
from Bio.PDB.Residue import Residue
from Bio.PDB.Vector import calc_dihedral

from Bio.PDB.TorusDBN import TorusDBNModel
from Bio.PDB.TorusDBN.TorusDBNExceptions import TorusDBNBuildPolypeptideException, \
    TorusDBNChainBreakException
from Bio.PDB.TorusDBN._utils import dssp_to_index

import numpy
import math
import time
import os



class TorusDBNTrainer(object):
    """
    Allow training a TorusDBN model with a given training set.
    """

    def __init__(self, seed=int(time.time()), show_info=False, show_warnings=True):
        # Mocapy config
        mocapy_seed(seed)
        
        # Print information and warnings to the standard output
        self.show_info = show_info
        self.show_warnings = show_warnings
                            
        # Inference engine parameters
        self.em_steps = 100
        self.burnin_steps = 20
        self.mcmc_steps = 1
        
        # Convergence check for the optimal node size search
        self.check_convergence = False
        self.convergence_window = 80
        self.convergence_threshold = 55

        # Default Model
        self.model = TorusDBNModel(seed)
        
        # Training sequences and mismasks 
        self.seq_list = []
        self.mismask_list = []
        self.__missing_residues = {}
            
            
    def get_model(self):
        return self.model
        
        
    def __read_missing_residues(self, filename):  
        
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
                            self.__missing_residues[
                                (chain_pdb, int(chain_index), int(res_index))
                            ] = Residue(id=('', int(res_index), ''), resname= res_name, segid='')
                    except ValueError:
                        self.info(
                            "Could not read missing residues file %s at line %d. Skipped." % 
                            (filename, i + 1)
                        )
        except IOError:
            self.info("Could not open file %s for reading missing residues." % filename)
        
        
    def train(self, training_set, missing_residues=None, use_aic=False):
        if missing_residues is not None:
            self.__read_missing_residues(missing_residues)
            
        self.seq_list, self.mismask_list = self.__create_sequence_and_mismask(
            training_set)   
        self.info('Training started...')
        ic = self.__train(use_aic)
        self.info('Training finished.')
        return ic        
        
        
    def __train(self, use_aic):        
        dbn = self.model.dbn
        
        self.hmm_ll_calculator = LikelihoodInfEngineHMM(
            dbn=dbn, hidden_node_index=0, check_dbn=False)
                    
        mcmc = GibbsRandom(dbn)
        em = EMEngine(dbn, mcmc, self.seq_list, self.mismask_list, [])
        
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
                if i > self.convergence_window:
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
            dbn=self.model.dbn, hidden_node_index=0, check_dbn=False)
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
            dbn=self.model.dbn, hidden_node_index=0, check_dbn=False)
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
        size_h = self.model.size_h
        return (size_h - 1) + size_h * (
            (size_h - 1) + parameters_d + (self.model.size_aa - 1) + 
            (self.model.size_ss - 1) + (self.model.size_cis - 1)
        )     
        
                
    def __create_sequence_and_mismask(self, training_set):
        seq_list = []
        mismask_list = []
        for filename in training_set:
            self.info('Reading data from training file %s...' % (filename))
            try:
                sequences, mismasks = self.__get_data_from_file(filename)
                seq_list += sequences
                mismask_list += mismasks
            except TorusDBNBuildPolypeptideException:
                self.info("Could not create polypeptide list from file %s ." % (filename))
            except TorusDBNChainBreakException, error:
                self.info(error.message)
                self.info("The file %s has chain breaks and was not included in the training set." % (filename))

        return seq_list, mismask_list
        
        
    def __get_missing_residue(self, chain_pdb, chain_index, residue_index):        
        chain_pdb = os.path.split(chain_pdb)[-1]
        try:
            return self.__missing_residues[(chain_pdb, chain_index, residue_index)]
        except KeyError:
            raise TorusDBNChainBreakException(
                "Chain break in file %s, chain %d, residue %d, could not be handled." % 
                (chain_pdb, chain_index, residue_index)
            )
        
        
    def __get_pp_with_chain_break(self, chain_pdb, pp_list):
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
                    pp_chains[-1].append(self.__get_missing_residue(
                        chain_pdb, len(pp_chains), res_index))
                pp_chains[-1] += pp

        return pp_chains, missing_residues


    def __get_conformation(self, res, prev_res):
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
        
        
    def __get_data_from_file(self, chain_pdb):
        sequences = []
        mismasks = []
                
        output_data = []
        output_mismask = []
        
        parser = PDBParser(QUIET=not self.show_warnings)
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
           pp_chains, chain_missing_residues = self.__get_pp_with_chain_break(chain_pdb, pp_list)
                       
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
                        seq[5] = self.__get_conformation(res, pp[i-1])

                output_data.append(seq)
                output_mismask.append(mism)

        sequences.append(numpy.array(output_data))
        mismasks.append(numpy.array(output_mismask, dtype = numpy.uint))
        return sequences, mismasks       
            
            
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
                    self.info("Training with node size = %d (sample %d)" % (i, j))
                    self.model.create_dbn(hidden_node_size=i)
                    IC = self.__train(use_aic)
                    IC_array[j].append(IC)
                    IC_cum += IC
                    
                    if (full_ll_dec):
                        # Save forward likelihoods in order to infer if it is decreasing
                        hmm_ll_calculator = LikelihoodInfEngineHMM(
                            dbn=self.model.dbn, hidden_node_index=0, check_dbn=False)
                        ll_full = hmm_ll_calculator.calc_ll(self.seq_list, self.mismask_list)
                        avg_full_LL[-1] = avg_full_LL[-1] + ll_full/self.__get_observation_count()
                    
                # Calculate mean IC for each node-size and add to array
                IC_array[node_samples + 1].append(IC_cum / node_samples)
                
                # Check if log-likelihood is decreasing 
                if (len(avg_full_LL) > 1) and (avg_full_LL[-1] < avg_full_LL[-2]) and \
                    (start_res == start_int) and full_ll_dec:
                    self.info("Log-likelihood is decreasing. There is no reason to test higher node sizes.")
                    break
                                     
            # Column number for maximum IC value
            max_position = IC_array[node_samples + 1].index(max(IC_array[node_samples + 1])) 
            self.info("Optimal node size: %s\n" % (IC_array[0][max_position]))
          
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
            self.model.create_dbn(hidden_node_size=IC_max_node)
            IC = self.__train(use_aic)
            IC_list.append(IC)
            dbn_list.append(self.model.dbn)
            
        IC_max = max(IC_list)
        self.model.dbn = dbn_list[IC_list.index(IC_max)]
        
        self.info("Optimal Model:\nHidden node size = %s\nIC = %s\n" % (IC_max_node, IC_max))        
        return IC_max_node, IC_max
        

    def info(self, message):
        if self.show_info:
            print(message)
