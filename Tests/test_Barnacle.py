import unittest
import os

from Bio.PDB.Barnacle.Barnacle import Barnacle


class BarnacleTestCase(unittest.TestCase):

    def setUp(self):
        current_directory = os.path.dirname(globals()["__file__"])
        self.datadir = os.path.join(current_directory, "Barnacle/")


    def test_barnacle_sample(self):
        # The target sequence
        sequence = "ACGU"

        # Initialize the model
        model = Barnacle(sequence, seed=123)

        # Draw first sample
        model.sample()

        # Sample position j
        j = 0       
        model.sample(start=j, end=j+1)
        
        obtained_filename = "test_Barnacle.pdb.obtained"
        expected_filename = "test_Barnacle.pdb.expected"

        # Save the structure
        model.save_structure(os.path.join(self.datadir, obtained_filename))

        # Compare generated file
        obtained_file = open(os.path.join(self.datadir, obtained_filename))
        expected_file = open(os.path.join(self.datadir, expected_filename))

        self.assertEquals(obtained_file.readlines(), expected_file.readlines())

    
    def test_barnacle_log_likelihood(self):
        # The target sequence
        sequence = "ACGU"

        # Initialize the model
        model = Barnacle(sequence, seed=123)

        # Draw first sample
        model.sample()

        # Sample position j
        j = 0       
        model.sample(start=j, end=j+1)

        # Calculate and print likelihood
        ll = model.get_log_likelihood()
        # print "test_Barnacle.pdb: ll = %f" % (ll)
        self.assertAlmostEquals(ll, 19.648623, places=4)
        
        

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
