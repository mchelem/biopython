# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os

from Bio.PDB import *
from Bio.PDB.TorusDBN.TorusDBN import TorusDBN


class TorusTestCase(unittest.TestCase):

    def get_training_set(self):
        return ['1BJQ.pdb', '1FAT.pdb'] #, '1MBN.pdb', '3QKS.pdb', '4MBN.pdb', '6GCH.pdb', '7DFR.pdb']
        
        
    def test_model_training(self):        
        training_set = self.get_training_set()        

        model = TorusDBN()
        
        BIC = model.train(training_set)
        self.assertAlmostEquals(BIC, -22374.2444164, places=4)
        
        AIC = model.train(training_set, use_aic=True)
        self.assertAlmostEquals(AIC, -10177.72654, places=4)

                
    def test_find_optimal_model(self):
        training_set = self.get_training_set()        
        
        model = TorusDBN()
        hidden_node_size, IC = model.find_optimal_model(training_set,
            node_samples=2, full_ll_dec=True)
            
        self.assertEquals(hidden_node_size, 5)
        self.assertAlmostEquals(IC , -1942.49842302, places=4)
            
          

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
