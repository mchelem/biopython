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
    
    
    def check_angles(self, angles_obtained, angles_expected):
        self.assertEquals(len(angles_obtained), len(angles_expected))
        
        for i in xrange(len(angles_obtained)):
            self.assertEquals(len(angles_obtained[i]), 2)
            self.assertAlmostEquals(angles_obtained[i][0], angles_expected[i][0])
            self.assertAlmostEquals(angles_obtained[i][1], angles_expected[i][1])
        
        
    def test_model_sample(self):       
        model = TorusDBN(seed=123)
        
        # Setting only the aa sequence
        model.set_aa('ACDEFGHIK')                     
        sample = model.sample()   
        
        self.assertEquals(model.get_aa(), 'ACDEFGHIK')
        self.assertEquals(model.get_ss(), 'ECHCECEEE')
        self.assertEquals(
            model.get_cis().tolist(), 
            [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0],
        )
        angles_expected = [
            [2.166881003185333, -1.6556858216141894], 
            [0.7400624722095404, 2.2578738132852965], 
            [-0.05282326657352576, -0.6030190954028775], 
            [1.938464071387846, 1.4833623359658787], 
            [-1.9887756273827366, -0.3490076143190232], 
            [-0.35597989981895195, -1.6498515994372953], 
            [-0.477415314483151, -1.951192331990111], 
            [2.25689353850571, 1.9343839253436967], 
            [-1.048644559246696, 0.4615353850695407]
        ]        
        self.check_angles(model.get_angles().tolist(), angles_expected)
             
        self.assertAlmostEquals(model.get_log_likelihood(), -73.07113, places=4)
        
        # Setting aa sequence, secondary structure and cis
        model = TorusDBN(seed=123)        
        model.set_aa('ACDE')
        model.set_cis([1, 1, 0, 0])              
        model.set_ss('ECCE')
        sample = model.sample()    
        
        self.assertEquals(model.get_aa(), 'ACDE')
        self.assertEquals(model.get_ss(), 'ECCE')
        self.assertEquals(model.get_cis().tolist(), [1.0, 1.0, 0.0, 0.0])  
        angles_expected = [
            [0.8482119396488295, -2.031278691606202], 
            [1.1989388179360942, -1.6305194474895712],
            [-0.9015472823190614, 3.055749238578653], 
            [-2.9424056925158375, 1.0225598309477402],
        ]
        self.check_angles(model.get_angles().tolist(), angles_expected)        
        
        self.assertAlmostEquals(model.get_log_likelihood(), -32.47606, places=4)
                
        # Setting aa sequence and angles
        model = TorusDBN(seed=123)   
        model.set_aa('ACDE')
        model.set_angles([[1.49, -0.64, 1.32, -0.66], [-1.64, 1.34, -0.79, 1.29]])        
        sample = model.sample()   
        
        self.assertEquals(model.get_aa(), 'ACDE')
        self.assertEquals(model.get_ss(), 'CCCH')
        self.assertEquals(model.get_cis().tolist(), [1.0, 0.0, 1.0, 0.0])  
        angles_expected = [
            [1.49, -1.64], 
            [-0.64, 1.34], 
            [1.32, -0.79], 
            [-0.66, 1.29],
        ]
        
        self.check_angles(model.get_angles().tolist(), angles_expected)                   
        self.assertAlmostEquals(model.get_log_likelihood(), -31.80456, places=4)        
        model.save_structure(os.path.join(self.datadir, 'test_TorusDBN.pdb'))                    
            


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
