# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

""" Methods to calculate coordinates for atoms given a set of angles. """

from math import sin, cos, pi, sqrt
import numpy

from Bio.PDB.Vector import Vector, rotaxis2m, calc_angle
from Bio.PDB.Polypeptide import index_to_three

from Bio.PDB.TorusDBN.TorusDBNExceptions import TorusDBNException
from Bio.PDB.TorusDBN._constants import ANGLE_N_CA_CB, ANGLE_CB_CA_C, \
    BOND_LENGTH_CA_CB, ANGLE_O_C_N, BOND_LENGTH_C_N, BOND_LENGTH_C_O, ATOMS, \
    get_bond_length, get_bond_angle
    
    

# Coordinate identifiers to index the coordinates matrix
X = 0
Y = 1
Z = 2    

# The number of atoms per residue, including CB and O
atoms_per_residue = 5 

# The number of atoms in the backbone
backbone_atoms = 3

  
def calculate_coordinates(angles, cis, aa_seq):
    """ Calculate coordinates given an amino acid residues, angles and cis/trans
    conformation.
    
    @param angles: Angles.
    @type angles: list(float)
    
    @param cis: Conformation of the peptide bonds.
    @type cis: list(int)
    
    @param aa_seq: Amino acid sequence.
    @type aa_seq: list(int)
    
    @rtype: numpy.array
    @return: The calculated coordinates.
    
    """
    coords = numpy.zeros((3, len(aa_seq) * atoms_per_residue))
        
    _place_backbone_atoms(coords, angles, cis, aa_seq)
    _place_CB_O_atoms(coords, angles, cis, aa_seq)
    
    return numpy.reshape(coords, coords.size,order='F').reshape(
        (len(aa_seq) * atoms_per_residue, 3))[:-1]
        
     
def _place_backbone_atoms(coords, angles, cis, aa_seq):
    """ Place atoms that are in the backbone.
    
    @param coords: Array where the coordinates are to be saved.
    @type coords: numpy.array
    
    @param angles: Angles.
    @type angles: list(float)
    
    @param cis: Conformation of the peptide bonds.
    @type cis: list(int)
    
    @param aa_seq: Amino acid sequence.
    @type aa_seq: list(int)
    
    """    
    for i in xrange(len(aa_seq)):
        current_res = index_to_three(aa_seq[i])
        current_atom = 'N'
        prev_atom =  'C'
        prev_prev_atom = 'CA'
        
        for j in xrange(backbone_atoms):
            index = i * atoms_per_residue + j
            if i == 0 and current_atom == 'N':
                coords[:,index] = 0
            elif i == 0 and current_atom == 'CA':
                coords[:,index] = 0
                coords[X][index] = get_bond_length(
                    current_atom, prev_atom, current_res)   
            else:
                v1 = _create_vector(coords, index, j, 1) 
                v2 = _create_vector(coords, index, j, 2) 
            
                if i == 0 and current_atom == 'C':
                    v3 = Vector(0, -1.0, 0)
                else:
                    v3 = _create_vector(coords, index, j, 3) 
                
                angle = get_bond_angle(
                    prev_prev_atom, prev_atom, current_atom, current_res)                
                dihedral = _calculate_dihedral(prev_atom, angles, cis, i)                
                bond_length = get_bond_length(prev_atom, current_atom, current_res)                 
                bond_length_prev = get_bond_length(
                    prev_prev_atom, prev_atom, current_res)                                    
                    
                D = Vector(
                    bond_length * cos(pi - angle),
                    bond_length * cos(pi - dihedral) * sin(pi - angle),
                    bond_length * sin(pi - dihedral) * sin(pi - angle),
                )                                                 
                bc = (v1 - v2) / bond_length_prev
                n = ((v1 - v3) ** bc).normalized()
                nbc = bc ** n
                basis_change = numpy.array((bc, nbc, n)).transpose()
                D = D.left_multiply(basis_change) + v1                
                coords[:,index] = D
                
            prev_prev_atom = prev_atom
            prev_atom = current_atom
            current_atom = ATOMS[(ATOMS.index(current_atom) + 1) % backbone_atoms]
                
    
def _place_CB_O_atoms(coords, angles, cis, aa_seq):
    """ Place CB and O atoms.
    
    @param coords: Array where the coordinates are to be saved.
    @type coords: numpy.array
    
    @param angles: Angles.
    @type angles: list(float)
    
    @param cis: Conformation of the peptide bonds.
    @type cis: list(int)
    
    @param aa_seq: Amino acid sequence.
    @type aa_seq: list(int)
    
    """
    L = len(aa_seq)
    for i in xrange(L):
        index = i * atoms_per_residue
        current_residue = index_to_three(aa_seq[i])
        
        # Place CB
        coords[:,index + 3] = _calculate_CB_coordinates(
            Vector(coords[:,index]),
            Vector(coords[:,index + 1]),
            Vector(coords[:,index + 2]),        
            current_residue,
        )
        
        # Place O
        if i < (L - 1):
            coords[:, index + 4] = _calculate_O_coordinates(
                Vector(coords[:,index + 1]),
                Vector(coords[:,index + 2]),  
                Vector(coords[:, index + atoms_per_residue]),
                current_residue,       
            )
            
        
def _calculate_CB_coordinates(position_N, position_CA, position_C, residue):
    """ Calculate coordinates for the CB atom.
    
    @param position_N: Coordinates of the N atom
    @type position_N: Vector
    
    @param position_CA: Coordinates of the CA atom
    @type position_CA: Vector
    
    @param position_C: Coordinates of the C atom
    @type position_C: Vector
    
    @param residue: Amino acid.
    @type residue: str
    
    @rtype: Vector
    @return: The calculated CB position
        
    """
    angle_N_CA_CB = ANGLE_N_CA_CB[residue] * pi / 180
    angle_CB_CA_C = ANGLE_CB_CA_C[residue] * pi / 180
    bond_length_CA_CB = BOND_LENGTH_CA_CB[residue]
        
    position_CB = _calculate_position_from_angles(
        position_N, position_CA, position_C, angle_N_CA_CB, 
        angle_CB_CA_C, bond_length_CA_CB
    )    
    if position_CB is None:
        CAN = position_N - position_CA
        CAC = position_C - position_CA
        
        # find rotation matrix that rotates n -120 degrees along the ca-c vector
        rotation_matrix = rotaxis2m(-pi * 120 / 180, CAC)
        CACB = CAN.left_multiply(rotation_matrix).normalized() * bond_length_CA_CB
        position_CB = CACB + position_CA
        
    return position_CB
    
    
def _calculate_O_coordinates(position_CA, position_C, position_N, residue):
    """ Calculate coordinates for the O atom.
    
    @param position_CA: Coordinates of the CA atom
    @type position_CA: Vector
    
    @param position_C: Coordinates of the C atom
    @type position_C: Vector
    
    @param position_N: Coordinates of the N atom
    @type position_N: Vector
    
    @param residue: Amino acid.
    @type residue: str
    
    @rtype: Vector
    @return: The calculated O position
        
    """
    angle_O_C_N = ANGLE_O_C_N[residue] * pi / 180
    bond_length_C_N = BOND_LENGTH_C_N[residue]
    
    D = Vector(
        BOND_LENGTH_C_O * cos(angle_O_C_N), BOND_LENGTH_C_O * sin(angle_O_C_N), 0)
    bc = (position_N - position_C) / bond_length_C_N
    n = ((position_C - position_CA) ** bc).normalized()
    nbc = bc ** n
    basis_change = numpy.array((bc, nbc, n)).transpose()
        
    return D.left_multiply(basis_change) + position_C
    
            
def _calculate_dihedral(atom, angles, cis, position):
    """ Calculate dihedral given atom, angles and current position.  
    
    @type atom: str
    @type angles: list(float)
    @type cis: int
    @type position: int
    """
    if atom == 'N':
        dihedral = pi if cis[position] == 0 else 0
    elif atom == 'CA':
        dihedral = angles[position][0]
    elif atom == 'C':
        dihedral = angles[position - 1][1]
    else:
        TorusDBNException("Could not calculate dihedral for atom %s." % (atom))
    
    return dihedral
       
    
def _calculate_position_from_angles(B, A, C, alpha, beta, l):
    """ Calculate position from three vector and two angles. 
    
        - B, A and C are the vectors    
        - alpha and beta are the angles   
        - l is the bond length.    
        
    @type A: Vector
    @type B: Vector
    @type C: Vector
    
    @type alpha: float
    @type beta: float
    @type l: float
    
    """
    
    gamma = calc_angle(B, A, C)
    Dx = l * cos(alpha)
    Dy = l * (cos(beta) - cos(alpha) * cos(gamma)) / sin(gamma)

    try:
        Dz = sqrt(l * l - Dx * Dx - Dy * Dy)
    except ValueError:
        return None
    
    Dab = Dx - cos(gamma) / sin(gamma) * Dy
    Dac = Dy / sin(gamma)
    Dabac = Dz
    
    AB = (B - A).normalized()
    AC = (C - A).normalized()
    ABxAC = (AB ** AC).normalized()  
    AD = AB * Dab + AC * Dac + ABxAC * Dabac
        
    return A + AD
    

def _create_vector(array, index, j, position):
    """ Method to create vectors for the backbone placement. """
    if (j - position) >= 0:
        array_index = index - position
    else: 
        array_index = index - atoms_per_residue + (backbone_atoms - position)
    
    return Vector(array[:, array_index])

