# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os


from Bio.PDB.TorusDBN.TorusDBNTrainer import TorusDBNTrainer



class TorusDBNTrainerTestCase(unittest.TestCase):
    
    def setUp(self):
        current_directory = os.path.dirname(globals()["__file__"])
        self.datadir = os.path.join(current_directory, "TorusDBN/")
    
    
    def get_training_set(self):
        pdb_dir = os.path.join(self.datadir, "PDB")
        pdb_files =  [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir)]
        return pdb_files
    
    
    def test_model(self):        
        training_set = self.get_training_set()  

        # setting seed for reproducibility
        trainer = TorusDBNTrainer(seed=123, show_info=True, show_warnings=False) 
        trainer.em_steps = 10 
        trainer.burnin_steps = 5
        
        BIC = trainer.train(
            training_set, 
            missing_residues=os.path.join(self.datadir, "missing_residues")
        )
        self.assertAlmostEquals(BIC, -63897.94317, places=4)
        
        AIC = trainer.train(training_set, use_aic=True)
        self.assertAlmostEquals(AIC, -36759.810827, places=4)   
        
        
    def test_model_optimization(self):
        training_set = self.get_training_set()        
        
        trainer = TorusDBNTrainer(seed=123, show_info=True, show_warnings=False)
        trainer.em_steps = 10
        trainer.burnin_steps = 5
        
        hidden_node_size, IC = trainer.find_optimal_model(training_set,
            node_samples=1, max_node=30, full_ll_dec=True)
            
        self.assertEquals(hidden_node_size, 5)
        self.assertAlmostEquals(IC ,  -1702.49251, places=4)
        
  
  
if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
