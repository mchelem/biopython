#  Barnacle: A probabilistic model of RNA conformational space
#
#  Copyright (C) 2008 Jes Frellsen, Ida Moltke and Martin Thiim 
#
#  Barnacle is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Barnacle is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Barnacle.  If not, see <http://www.gnu.org/licenses/>.

from mocapy.framework import DBN, NodeFactory, mocapy_seed
from mocapy.vonmises import VonMisesDensities
import model_parameters as param
import time


# Set constant values
num_angles = 7
h_size = 20

angle_id_pos = 0
hd_pos = 1
angle_pos = 2
num_nodes = 3


def make_model():
    """
    Method for constructing the model
    """    
    # Angle identifier
    angleId0 = NodeFactory.new_discrete_node(
        node_size=num_angles, name='angleId0', cpd=param.angleId0_cpd)
    angleId1 = NodeFactory.new_discrete_node(
        node_size=num_angles, name='angleId1', cpd=param.angleId1_cpd)

    # Hidden node
    hd0 = NodeFactory.new_discrete_node(
        node_size=h_size, name='hd0', cpd=param.hd0_cpd)
    hd1 = NodeFactory.new_discrete_node(
        node_size=h_size, name='hd1', cpd=param.hd1_cpd)

    # Angle node
    mus = param.angle_mus
    mus.shape = 20
    mus = mus.tolist()
    kappas = param.angle_kappas.tolist()

    angle = NodeFactory.new_vonmises_node(
        name='angle', 
        mus=mus,
        kappas=kappas, 
    )

    start_nodes = [angleId0, hd0, angle]
    end_nodes = [angleId1, hd1, angle]
    dbn = DBN(start_nodes, end_nodes)

    # Internal connections
    dbn.add_intra('angleId0', 'hd0')
    dbn.add_intra('hd0','angle')

    # Connections between consecutive slices
    dbn.add_inter('hd0','hd1')
    dbn.add_inter('angleId0','angleId1')

    dbn.construct()

    return dbn
