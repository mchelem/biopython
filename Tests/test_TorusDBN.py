# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os

from Bio.PDB import *
from Bio.PDB.TorusDBN.TorusDBN import TorusDBN


class BarnacleTestCase(unittest.TestCase):

    def get_training_set(self):
        return ['pdb1fat.ent']
        

    def test_model_sample(self):        
        training_set = self.get_training_set()        
        
        model = TorusDBN()
        model.train(training_set)
        
        print 'AIC', model.calculate_AIC()
        print 'BIC', model.calculate_BIC()
        
        

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
