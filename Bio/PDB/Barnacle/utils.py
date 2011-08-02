# Copyright (C) 2008 by Jes Frellsen, Ida Moltke and Martin Thiim
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import numpy
from numpy import array, zeros, ones
from model import angle_id_pos, angle_pos, num_nodes, num_angles


########## FUNCTION FOR MAKING DATA AND MISMASKS ##########
angle_ids = [0.0, 1.0, 2.0, 6.0, 3.0, 4.0, 5.0]

def make_data_and_mism(len_nucs):
    """
    Make a sequence and a mismask for model of length 'len_nucs'
    """
    len_angles = len_nucs*num_angles

    data = zeros((len_angles, num_nodes))
    mism_sample = ones((len_angles, num_nodes), dtype=numpy.uint)
    mism = ones((len_angles, num_nodes), dtype=numpy.uint)

    # Set the angle IDs
    data[:,angle_id_pos] = array(angle_ids*len_nucs)
    mism_sample[:,angle_id_pos] = 0
    
    mism[:,angle_id_pos] = 0
    mism[:,angle_pos] = 0

    return data, mism_sample, mism
    

########## FUNCTION FOR CONVERTING DATA INTO AN ANGLELIST ##########

def data_to_anglelist(data):
    length = data.shape[0]
    angle_list = map(lambda l: tuple(map(float, list(l))), data[:,angle_pos].reshape((length/7,7))[:,(0,1,2,4,5,6,3)])
    return angle_list
