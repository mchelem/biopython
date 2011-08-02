# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os

from Bio.PDB import *
from Bio.PDB.TorusDBN.TorusDBN import TorusDBN


class TorusTestCase(unittest.TestCase):

    
    def setUp(self):
        current_directory = os.path.dirname(globals()["__file__"])
        self.datadir = os.path.join(current_directory, "TorusDBN/")
    
    
    def get_training_set(self):
        pdb_dir = os.path.join(self.datadir, "PDB")
        pdb_files =  [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir)]
        return pdb_files
    
    
    def test_find_optimal_model(self):
        training_set = self.get_training_set()        
        
        model = TorusDBN(seed=123)
        hidden_node_size, IC = model.find_optimal_model(training_set,
            node_samples=2, max_node=30, full_ll_dec=True)
            
        self.assertEquals(hidden_node_size, 5)
        self.assertAlmostEquals(IC , -1892.01090, places=4)
        
        
    def test_model_training(self):        
        training_set = self.get_training_set()      
        model = TorusDBN(seed=123)  # setting seed for reproducibility
        
        BIC = model.train(training_set)
        self.assertAlmostEquals(BIC, -22375.04932, places=4)
        
        AIC = model.train(training_set, use_aic=True)
        self.assertAlmostEquals(AIC, -10210.10216, places=4)
        
        
    def test_model_sample(self):       
        model = TorusDBN(seed=123)
        
        model.set_aa('ACDEFGHIK')        
        sample = model.sample()        
        self.assertAlmostEquals(model.get_log_likelihood(), -74.50195, places=4)
        
        model = TorusDBN(seed=123)        
        model.set_aa('ACDE')
        model.set_cis([1, 1, 0, 0])              
        model.set_ss([1, 2, 2, 1])
        sample = model.sample()        
        self.assertAlmostEquals(model.get_log_likelihood(), -33.111980, places=4)
                
        model = TorusDBN(seed=123)   
        model.set_aa('ACDE')
        model.set_angles([[1.49, -0.64, 1.32, -0.66], [-1.64, 1.34, -0.79, 1.29]])        
        sample = model.sample()        
        self.assertAlmostEquals(model.get_log_likelihood(), -27.74002, places=4)
            
          

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
