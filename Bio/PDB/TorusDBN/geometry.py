# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from Bio.PDB.Vector import Vector, rotaxis2m, calc_angle
from Bio.PDB.Polypeptide import index_to_three

from Bio.PDB.TorusDBN.EnghHuberConstants import ANGLE_N_CA_CB, ANGLE_CB_CA_C, \
    BOND_LENGTH_CA_CB, ANGLE_O_C_N, BOND_LENGTH_C_N, BOND_LENGTH_C_O, ATOMS
from Bio.PDB.TorusDBN.EnghHuberConstants import get_bond_length, \
    get_bond_angle
    
import math
import numpy

numpy.set_printoptions(precision=7, suppress=True)

X = 0
Y = 1
Z = 2    

atoms_per_residue = 5 # including CB and O
backbone_atoms = 3

  
def calculate_dihedral(prev_atom, angles, cis, i):
    if prev_atom == 'N':
        if cis[i] == 0:
            dihedral = math.pi
        else:
            dihedral = 0
    elif prev_atom == 'CA':
        dihedral = angles[i][0]
    elif prev_atom == 'C':
        dihedral = angles[i-1][1]
    else:
        ValueError("Illegal atom.")
    
    return dihedral
    

def create_vector(array, index, j, position):
    if (j - position) >= 0:
        array_index = index - position
    else: 
        array_index = index - atoms_per_residue + (backbone_atoms - position)
    
    v = Vector(
        array[X][array_index], 
        array[Y][array_index], 
        array[Z][array_index],
    )
    return v
    
    
def calculate_position_from_angles(B, A, C, alpha, beta, l, choice=0):
    
    gamma = calc_angle(B, A, C)
    Dx = l * math.cos(alpha)
    Dy = l * (math.cos(beta) - math.cos(alpha) * math.cos(gamma)) / math.sin(gamma)

    try:
        Dz = math.sqrt(l * l - Dx * Dx - Dy * Dy)
    except ValueError:
        return None
    
    Dab = Dx - math.cos(gamma) / math.sin(gamma) * Dy
    Dac = Dy / math.sin(gamma)
    Dabac = Dz
    
    AB = (B - A).normalized()
    AC = (C - A).normalized()
    ABxAC = (AB ** AC).normalized()
    
    if choice == 0:
        AD = AB * Dab + AC * Dac + ABxAC * Dabac
    else:
        AD = AB * Dab + AC * Dac - ABxAC * Dabac
        
    return A + AD


def place_CB_atom(position_N, position_CA, position_C, residue):
    angle_N_CA_CB = ANGLE_N_CA_CB[residue] * math.pi / 180
    angle_CB_CA_C = ANGLE_CB_CA_C[residue] * math.pi / 180
    bond_length_CA_CB = BOND_LENGTH_CA_CB[residue]
        
    position_CB = calculate_position_from_angles(
        position_N,
        position_CA, 
        position_C,
        angle_N_CA_CB,
        angle_CB_CA_C,
        bond_length_CA_CB
    )    
    if position_CB is None:
        CAN = position_N - position_CA
        CAC = position_C - position_CA
        
        # find rotation matrix that rotates n -120 degrees along the ca-c vector
        rotation_matrix = rotaxis2m(-math.pi * 120 / 180, CAC)
        CACB = CAN.left_multiply(rotation_matrix).normalized() * bond_length_CA_CB
        position_CB = CACB + position_CA
        
    return position_CB
    
    
def place_O_atom(position_CA, position_C, position_N, residue):
    angle_O_C_N = ANGLE_O_C_N[residue] * math.pi / 180
    bond_length_C_N = BOND_LENGTH_C_N[residue]
    
    D = Vector(
        BOND_LENGTH_C_O * math.cos(angle_O_C_N),
        BOND_LENGTH_C_O * math.sin(angle_O_C_N),
        0,
    )
    bc = (position_N - position_C) / bond_length_C_N
    n = ((position_C - position_CA) ** bc).normalized()
    nbc = bc ** n
    basis_change = numpy.array((bc, nbc, n)).transpose()
        
    return D.left_multiply(basis_change) + position_C
    
    
def get_coordinates_from_angles(angles, cis, aa):
    
    L = len(aa)
    
    coords = numpy.zeros((3, L * atoms_per_residue))
        
    for i in xrange(L):
        current_residue = index_to_three(aa[i])
        current_atom = 'N'
        prev_atom =  'C'
        prev_prev_atom = 'CA'
        
        for j in xrange(backbone_atoms):
            index = i * atoms_per_residue + j
            if i == 0 and current_atom == 'N':
                coords[X][index] = 0
                coords[Y][index] = 0
                coords[Z][index] = 0
            elif i == 0 and current_atom == 'CA':
                coords[X][index] = get_bond_length(
                    current_atom, prev_atom, current_residue)
                coords[Y][index] = 0
                coords[Z][index] = 0      
            else:
                v1 = create_vector(coords, index, j, 1) 
                v2 = create_vector(coords, index, j, 2) 
            
                if i == 0 and current_atom == 'C':
                    v3 = Vector(0, -1.0, 0)
                else:
                    v3 = create_vector(coords, index, j, 3) 
                
                angle = get_bond_angle(
                    prev_prev_atom, prev_atom, current_atom, current_residue)
                
                dihedral = calculate_dihedral(prev_atom, angles, cis, i)
                
                bond_length = get_bond_length(
                    prev_atom, current_atom, current_residue)
                    
                
                    
                D = Vector(
                    bond_length * math.cos(math.pi - angle),
                    bond_length * math.cos(math.pi - dihedral) * math.sin(math.pi - angle),
                    bond_length * math.sin(math.pi - dihedral) * math.sin(math.pi - angle),
                )
                
                bond_length_prev = get_bond_length(
                    prev_prev_atom, prev_atom, current_residue)
                                    
                bc = (v1 - v2) / bond_length_prev
                n = ((v1 - v3) ** bc).normalized()
                nbc = bc ** n
                basis_change = numpy.array((bc, nbc, n)).transpose()
                D = D.left_multiply(basis_change) + v1
                
                coords[X][index] = D[X]
                coords[Y][index] = D[Y]
                coords[Z][index] = D[Z]
                
            prev_prev_atom = prev_atom
            prev_atom = current_atom
            current_atom = ATOMS[(ATOMS.index(current_atom) + 1) % backbone_atoms]
                
    
    for i in xrange(L):
        index = i * atoms_per_residue
        current_residue = index_to_three(aa[i])
        
        # Place CB
        position_CB = place_CB_atom(
            Vector(coords[X][index], coords[Y][index], coords[Z][index]),
            Vector(coords[X][index + 1], coords[Y][index + 1], coords[Z][index + 1]),
            Vector(coords[X][index + 2], coords[Y][index + 2], coords[Z][index + 2]),        
            current_residue,
        )
        coords[X][index + 3] = position_CB[X]
        coords[Y][index + 3] = position_CB[Y]
        coords[Z][index + 3] = position_CB[Z]
        
        # Place O
        if i < (L - 1):
            position_O = place_O_atom(
                Vector(coords[X][index + 1], coords[Y][index + 1], coords[Z][index + 1]),
                Vector(coords[X][index + 2], coords[Y][index + 2], coords[Z][index + 2]),
                Vector(
                    coords[X][index + atoms_per_residue], 
                    coords[Y][index + atoms_per_residue], 
                    coords[Z][index + atoms_per_residue],        
                ),
                current_residue,       
            )
            coords[X][index + 4] = position_O[X]
            coords[Y][index + 4] = position_O[Y]
            coords[Z][index + 4] = position_O[Z]
            
        else:
            #FIXME: No O placement at last position
            pass
            
    return numpy.reshape(coords, coords.size,order='F').reshape(
        (L * atoms_per_residue, 3))[:-1]

