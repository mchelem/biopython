# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import numpy
import math
import os
import warnings
from time import time

from mocapy.framework import EMEngine, mocapy_seed, eMISMASK
from mocapy.inference import GibbsRandom, LikelihoodInfEngineHMM

from Bio.PDB.TorusDBN import TorusDBNModel
from Bio.PDB.TorusDBN.TorusDBNExceptions import TorusDBNBuildPolypeptideException, \
    TorusDBNChainBreakException, TorusDBNWarning
from Bio.PDB.TorusDBN._io import read_missing_residues, create_sequence_from_file



class TorusDBNTrainer(object):
    """ Allow training a TorusDBN model with a given training set. """
    

    def __init__(self, seed=int(time()), show_info=False, show_warnings=True):
        """
        @param seed: Seed for the random number generator.
        @type seed: int
        
        @param show_info: Output information during the training.
        @type show_info: bool
        
        @param show_warnings: Output warnings.
        @type show_warnings: bool
        
        """
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
                            
        
    def train(self, training_set, missing_residues=None, use_aic=False):
        """ Train the model with a given training set.
        
        @param training_set: The files to be used in the training.
        @type training_set: list(str)
        
        @param missing_residues: The file describing the missing residues in
            the training set.
        @type missing_residues: str
        
        @param use_aic: Use the Akaike information criterion to measure the 
            model fitness.
        @type use_aic: bool
        
        """
        if not self.show_warnings:
            warning_list = warnings.filters[:]
            warnings.filterwarnings('ignore', category=TorusDBNWarning)
            
        if missing_residues is not None:
            missing_residues = read_missing_residues(missing_residues)
            
        self.seq_list, self.mismask_list = self._create_sequence_and_mismask(
            training_set, missing_residues)         
        self.info('Training started...')
        ic = self._train(use_aic)
        self.info('Training finished.')
        
        if not self.show_warnings:
            warnings.filters = warning_list
            
        return ic        
        
        
    def _train(self, use_aic):    
        """ Train the model.
        
        @param use_aic: Use the Akaike information criterion to measure the 
            model fitness.
        @type use_aic: bool
        
        """    
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
                     if self._get_hairiness(ll_list):
                        break
                        
        if (use_aic):
            return self.calculate_AIC()
        else:
            return self.calculate_BIC()

    
    def calculate_AIC(self):
        """ The Akaike information criterion (AIC) is a measure of the relative 
        goodness of fit of a statistical model.
        
        Akaike, Hirotugu (1974). "A new look at the statistical model 
        identification". IEEE Transactions on Automatic Control 19 (6): 
        716-723. doi:10.1109/TAC.1974.1100705. MR0423716.
        """            
        hmm_ll_calculator = LikelihoodInfEngineHMM(
            dbn=self.model.dbn, hidden_node_index=0, check_dbn=False)
        ll_full = hmm_ll_calculator.calc_ll(self.seq_list, self.mismask_list)
        return 2 * ll_full - 2 * self._get_parameter_count()
        
                
    def calculate_BIC(self):
        """ The Bayesian information criterion (BIC) is a criterion for model 
        selection among a finite set of models.
        
        Schwarz, Gideon E. (1978). "Estimating the dimension of a model". Annals
        of Statistics 6 (2): 461-464. doi:10.1214/aos/1176344136. MR468014.
        """            
        hmm_ll_calculator = LikelihoodInfEngineHMM(
            dbn=self.model.dbn, hidden_node_index=0, check_dbn=False)
        ll_full = hmm_ll_calculator.calc_ll(self.seq_list, self.mismask_list)     
        return 2 * ll_full - self._get_parameter_count() * math.log(
            self._get_observation_count())

    
    def __get_hairiness(self, ll_list):
        """ Check whether likelihood is oscilating during the training.
        
        @param ll_list: Likelihood list.
        @type ll_list: list(float)
        
        """
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
                
                
    def _create_sequence_and_mismask(self, training_set, missing_residues):
        """ Create sequence and mismask (mask that identifies values as
            hidden, observed or missing) from a training set.
            
        @param training_set: The files to be used in the training.
        @type training_set: list(str)
        
        @param missing_residues: The file describing the missing residues in
            the training set.
        @type missing_residues: str
        
        """
        seq_list = []
        mismask_list = []
        
        if not self.show_warnings:
            warning_list = warnings.filters[:]
            warnings.filterwarnings('ignore', category=TorusDBNWarning)
            
        for filename in training_set:
            self.info('Reading data from training file %s...' % (filename))
            try:
                sequences, mismasks = create_sequence_from_file(
                    filename, missing_residues)
                seq_list += sequences
                mismask_list += mismasks
            except (TorusDBNBuildPolypeptideException, 
                TorusDBNChainBreakException) as error:
                warnings.warn(
                    "%s The file was not included in the training set." % error,
                    TorusDBNWarning
                )
                
        if not self.show_warnings:
            warnings.filters = warning_list
            
        return seq_list, mismask_list        
            
            
    def find_optimal_model(
        self, training_set, use_aic=False, min_node=10, 
        max_node=90, start_size=20, end_size=5, node_samples=4, 
        check_decreasing_ll=False, missing_residues=None):
        """
        Optimization method to find the best size for the hidden node 
        according to a training set.
        
        @param training_set: The files to be used in the training.
        @type training_set: list(str)
        
        @param min_node: Minimum size of the hidden node
        @type min_node: int
        
        @param max_node: Maximum size of the hidden node
        @type max_node: int
        
        @param start_size: Start size of the hidden node
        @type start: int
        
        @param end_size: Final size of the hidden node
        @type end: int
        
        @param use_aic: Use the Akaike information criterion to measure the 
            model fitness.
        @type use_aic: bool        
        
        @param check_decreasing_ll: Check whether the loglikelihood is decreasing.
        @type check_decreasing_ll: bool
        
        @param missing_residues: The file describing the missing residues in
            the training set.
        @type missing_residues: str
        
        
        """            
        if not self.show_warnings:
            warning_list = warnings.filters[:]
            warnings.filterwarnings('ignore', category=TorusDBNWarning)
            
        if missing_residues is not None:
            missing_residues = read_missing_residues(missing_residues)
            
        self.seq_list, self.mismask_list = self._create_sequence_and_mismask(
            training_set, missing_residues)             
            
        max_position = 0
        start_res = start_size
        avg_full_LL = []
        
        IC_array = [[]*n for n in xrange(node_samples + 2)]
        
        # Decrease size resolution until threshold (end_size)
        while start_size >= end_size:
            # Loop over node sizes
            for i in xrange(min_node, max_node + 1, start_size):
                
                # Continues if at the maximum node size from the previous resolution 
                if (len(IC_array[0]) > 0 and i == IC_array[0][max_position]) or i <= 0:
                    continue

                # Add node-size value to header
                IC_array[0].append(i)
                IC_cum = 0
                
                if start_res == start_size:
                    avg_full_LL.append(0)
                    
                for j in xrange(1, node_samples + 1):
                    self.info("Training with node size = %d (sample %d)" % (i, j))
                    self.model.create_dbn(hidden_node_size=i)
                    IC = self._train(use_aic)
                    IC_array[j].append(IC)
                    IC_cum += IC
                    
                    if (check_decreasing_ll):
                        # Save forward likelihoods in order to infer if it is decreasing
                        hmm_ll_calculator = LikelihoodInfEngineHMM(
                            dbn=self.model.dbn, hidden_node_index=0, check_dbn=False)
                        ll_full = hmm_ll_calculator.calc_ll(self.seq_list, self.mismask_list)
                        avg_full_LL[-1] = avg_full_LL[-1] + ll_full/self._get_observation_count()
                    
                # Calculate mean IC for each node-size and add to array
                IC_array[node_samples + 1].append(IC_cum / node_samples)
                
                # Check if log-likelihood is decreasing 
                if (len(avg_full_LL) > 1) and (avg_full_LL[-1] < avg_full_LL[-2]) and \
                    (start_res == start_size) and check_decreasing_ll:
                    self.info("Log-likelihood is decreasing. There is no reason to test higher node sizes.")
                    break
                                     
            # Column number for maximum IC value
            max_position = IC_array[node_samples + 1].index(max(IC_array[node_samples + 1])) 
            self.info("Optimal node size: %s\n" % (IC_array[0][max_position]))
          
            # Update resolution
            start_size = start_size / 2
            
            # Update node limits
            min_node = IC_array[0][max_position] - start_size
            max_node = IC_array[0][max_position] + start_size
     
        IC_max_node = IC_array[0][max_position]
        
        # Final train to the optimal model
        dbn_list = []
        IC_list = []
        
        for j in xrange(node_samples):
            self.model.create_dbn(hidden_node_size=IC_max_node)
            IC = self._train(use_aic)
            IC_list.append(IC)
            dbn_list.append(self.model.dbn)
            
        IC_max = max(IC_list)
        self.model.dbn = dbn_list[IC_list.index(IC_max)]
        
        self.info("Optimal Model:\nHidden node size = %s\nIC = %s\n" % (IC_max_node, IC_max)) 
        
        if not self.show_warnings:
            warnings.filters = warning_list       
        return IC_max_node, IC_max
        

    def info(self, message):
        """ Print a message to the standard output, in case show_info is enabled.
        
        @param message: The message to be printed.
        @type message: str
        
        """
        if self.show_info:
            print(message)
            
    
    def get_model(self):
        """ Get the current TorusDBN model. 
        
        @rtype: TorusDBNModel
        @return: The model.
                
        """
        return self.model
        
    
    def _get_observation_count(self):
        """ The total number of observations.
        
        @rtype: int
        @return: Number of observations.
        
        """
        observation_count = 0
        for sequence in self.seq_list:
            observation_count += sequence.shape[0]     
            
        return observation_count
        
                
    def _get_parameter_count(self):
        """ The number of parameters.
        
        @rtype: int
        @return: Number of parameters.
        
        """
        parameters_d = 5;
        size_h = self.model.size_h
        return (size_h - 1) + size_h * (
            (size_h - 1) + parameters_d + (self.model.size_aa - 1) + 
            (self.model.size_ss - 1) + (self.model.size_cis - 1)
        )     
