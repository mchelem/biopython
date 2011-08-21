# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os

from Bio.PDB.TorusDBN import TorusDBNModel



class TorusDBNModelTestCase(unittest.TestCase):
    
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
        model = TorusDBNModel(seed=123)
        
        # Setting only the aa sequence
        model.set_aa('ACDEFGHIK')                     
        model.sample()   
        
        self.assertEquals(model.get_aa(), 'ACDEFGHIK')
        self.assertEquals(model.get_ss(), 'ECCHHEHEE')
        self.assertEquals(
            model.get_cis().tolist(), 
            [1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0],
        )
        angles_expected = [
            [0.7946544230502923, -0.24202043516862584], 
            [0.1712591696301713, -1.6190874070629098], 
            [0.81024195153348, 0.17085127914513104], 
            [-0.09738088478327811, -1.643003405419636], 
            [-0.5527929198676731, -0.8342016939985792], 
            [0.5794731030215697, 1.3213757980610115], 
            [0.3150308170378142, 0.5547944828962164], 
            [-0.0725375713273676, -1.8327574409710987], 
            [-2.9424056925158375, 1.0225598309477402],
        ]
        self.check_angles(model.get_angles().tolist(), angles_expected)
             
        self.assertAlmostEquals(model.get_log_likelihood(), -73.07113, places=4)
        
        # Setting aa sequence, secondary structure and cis
        model = TorusDBNModel(seed=123)        
        model.set_aa('ACDE')
        model.set_cis([1, 1, 0, 0])              
        model.set_ss('ECCE')
        model.sample()    
        
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
        model = TorusDBNModel(seed=123)   
        model.set_aa('ACDE')
        model.set_angles([
            [1.49, -1.64],
            [-0.64, 1.34],
            [1.32, -0.79],
            [-0.66, 1.29],
        ])
        model.sample()   
    
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
            
    
    def test_save_load_dbn(self):
        filename = os.path.join(self.datadir, 'test_save.dbn')
        model = TorusDBNModel(seed=123)
        model.create_dbn(hidden_node_size=55)
        model.save_dbn(filename)        
        model.set_aa('ACDE')
        model.sample()

        loaded_model = TorusDBNModel(seed=123)
        loaded_model.load_dbn(filename)
        loaded_model.set_aa('ACDE')
        loaded_model.sample()
        
        self.check_angles(model.get_angles(), loaded_model.get_angles())        
        self.assertEquals(model.get_cis().tolist(), loaded_model.get_cis().tolist())
        self.assertEquals(model.get_ss(), loaded_model.get_ss())

        

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
