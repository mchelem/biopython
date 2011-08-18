# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

""" Bond length and angle constants.
These constants are primarily extracted from Engh, Huber, 1991.

Engh R A & Huber R (1991). Accurate bond and angle parameters for X-ray protein 
structure refinement. Acta Cryst., A47, 392-400. DOI: 10.1107/S0108767391001071
    
"""
import math



class DefaultValueDict(dict):
    """ Dictionary with default values for missing keys. """
    
    def __init__(self, default, values=None, **kwds):
        """
        @param default: The default value when the key is missing.
        @type: object
        
        @param values: The values to initialize the dictionary.
        @type: list or dict            
        """
        if values:
            super(DefaultValueDict, self).__init__(values)
        self.update(kwds)
        self.default = default


    def __getitem__(self, key):
        return self.get(key, self.default)


    def __copy__(self):
        return DefaultValueDict(self.default, self)


ATOMS = [
    'N', 'CA', 'C', 'O', 'CB', 'SG', 'OG', 'OG1', 'CG', 'CG1', 
    'CG2', 'SD', 'OD1', 'OD2', 'ND1', 'ND2', 'CD', 'CD1', 'CD2', 'OE1', 
    'OE2', 'NE', 'NE1', 'NE2', 'CE', 'CE1', 'CE2', 'CE3', 'NZ', 'CZ', 
    'CZ2', 'CZ3', 'OH', 'NH1', 'NH2', 'CH2', 'CM', 'H', 'HA', 'HA1', 
    'HA2', 'HA3', 'HB', 'HB1', 'HB2', 'HB3', 'HG', 'HG1', 'HG2', 'HG3', 
    'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'HD1', 'HD2', 'HD3', 
    'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23', 'HE', 'HE1', 'HE2', 
    'HE3', 'HE21', 'HE22', 'HH', 'HH2', 'HH11', 'HH12', 'HH21', 'HH22', 
    'HZ', 'HZ1', 'HZ2', 'HZ3', 'H1', 'H2', 'H3', 'OXT', 'XX', 'XS', 
    'XO', 'XN', 'XC', 'XH',
]
    
ANGLE_N_CA_CB = DefaultValueDict(110.5, {'ALA':110.4, 'ILE':111.5, 'THR':111.5, 'VAL':111.5,'PRO':103.0})
ANGLE_CB_CA_C = DefaultValueDict(110.1, {'ALA':110.5, 'ILE':109.1, 'THR':109.1, 'VAL':109.1})
BOND_LENGTH_CA_CB = DefaultValueDict(1.530, {'ALA':1.521, 'ILE':1.540, 'THR':1.540, 'VAL':1.540})
ANGLE_O_C_N = DefaultValueDict(123.0, {'PRO':122.0})
BOND_LENGTH_C_N = DefaultValueDict(1.329, {'PRO':1.341})
BOND_LENGTH_C_O = 1.231

BOND_LENGTH_MAINCHAIN = {
    ('CA', 'N'): DefaultValueDict(1.458, {'GLY': 1.451, 'PRO': 1.466}),
    ('C', 'CA'): DefaultValueDict(1.525, {'GLY': 1.516}),
    ('C', 'N'): DefaultValueDict(1.329, {'PRO': 1.341}),
    ('CA', 'CB'): DefaultValueDict(1.530, {'ALA': 1.521, 'ILE': 1.540, 'THR': 1.540, 'VAL': 1.540}),
    ('C', 'O'): DefaultValueDict(1.231),
    ('C', 'OXT'): DefaultValueDict(1.231),
    ('CA', 'CA'): DefaultValueDict(3.8),
}

BOND_LENGTH_SIDECHAIN = {
    'CYS': {('CB', 'SG'): 1.808},
    'ASP': {('CB', 'CG'): 1.516, ('CG', 'OD1'): 1.249, ('CG', 'OD2'): 1.249},
    'GLU': {('CB', 'CG'): 1.520, ('CD', 'CG'): 1.516, ('CD', 'OE1'): 1.249, 
        ('CD', 'OE2'): 1.249},
    'PHE': {('CB', 'CG'): 1.502, ('CD1', 'CG'): 1.384, ('CD2', 'CG'): 1.384, 
        ('CD1', 'CE1'): 1.382, ('CD2', 'CE2'): 1.382, ('CE1', 'CZ'): 1.382, 
        ('CE2', 'CZ'): 1.382},
    'HIS': {('CB', 'CD'): 1.497, ('CG', 'ND1'): 1.371, ('CD2', 'CG'): 1.356, 
        ('CE1', 'ND1'): 1.319, ('CD2', 'NE2'): 1.374, ('CE1', 'NE2'): 1.374},
    'ILE': {('CB', 'CG1'): 1.530, ('CB', 'CG2'): 1.521, ('CD1', 'CG1'): 1.513},
    'LYS': {('CB', 'CG'): 1.520, ('CD', 'CG'): 1.520, ('CD', 'CE'): 1.520, 
        ('CE', 'NZ'): 1.489},
    'LEU': {('CB', 'CG'): 1.530, ('CD1', 'CG'): 1.521, ('CD2', 'CG'): 1.521},
    'MET': {('CB', 'CG'): 1.520, ('CG', 'SD'): 1.803, ('CE', 'SD'): 1.791},
    'ASN': {('CB', 'CG'): 1.520, ('CG', 'OD1'): 1.231, ('CG', 'ND2'): 1.328},
    'PRO': {('CB', 'CG'): 1.492, ('CD', 'CG'): 1.503},
    'GLN': {('CB', 'CG'): 1.520, ('CD', 'CG'): 1.516, ('CD', 'OE1'): 1.231, 
        ('CD', 'NE2'): 1.328},
    'ARG': {('CB', 'CG'): 1.520, ('CD', 'CG'): 1.520, ('CD', 'NE'): 1.460, 
        ('CZ', 'NE'): 1.329, ('CZ', 'NH1'): 1.326, ('CZ', 'NH2'): 1.326},
    'SER': {('CB', 'OG'): 1.417}, 
    'THR': {('CB', 'OG1'): 1.433, ('CB', 'CG2'): 1.521},    
    'VAL': {('CB', 'CG1'): 1.521, ('CB', 'CG2'): 1.521},
    'TRP': {('CB', 'CG'): 1.498, ('CD1', 'CG'): 1.433, ('CD2', 'CG'): 1.365, 
        ('CD1', 'NE1'): 1.374, ('CE2', 'NE1'): 1.370, ('CD2', 'CE2'): 1.409, 
        ('CD2', 'CE3'): 1.398, ('CE2', 'CZ2'): 1.394, ('CH2', 'CZ2'): 1.368,
        ('CE3', 'CZ3'): 1.382, ('CH2', 'CZ3'): 1.400},
    'TYR': {('CB', 'CG'): 1.512, ('CD1', 'CG'): 1.389, ('CD2', 'CG'): 1.389, 
        ('CD1', 'CE1'): 1.382, ('CD2', 'CE2'): 1.382, ('CE1', 'CZ'): 1.378,
        ('CE2', 'CZ'): 1.378,  ('OH', 'CZ'): 1.376},  
}

BOND_LENGTH_PSEUDO_SIDECHAIN = {'ALA': 1.54, 'CYS': 2.8, 'ASP': 2.92, 'GLU': 3.125,
    'PHE': 3.79, 'GLY': 1.0, 'HIS': 3.57, 'ILE': 2.7, 'LYS': 4.6, 'LEU': 3.05,
    'MET': 3.185, 'ASN': 2.91, 'PRO': 2.29, 'GLN': 3.875, 'ARG': 4.8, 'SER': 2.43,
    'THR': 2.17, 'VAL': 2.19, 'TRP': 4.2, 'TYR': 4.27}

BOND_ANGLE_MAINCHAIN = {
    ('C', 'N', 'CA'): DefaultValueDict(121.7, {'GLY': 120.6, 'PRO': 122.6}),
    ('N', 'CA', 'C'): DefaultValueDict(111.2, {'GLY': 112.5, 'PRO': 111.8}),
    ('CA', 'C', 'N'): DefaultValueDict(116.2, {'GLY': 116.4, 'PRO': 1116.9}),
    ('CA', 'C', 'O'): DefaultValueDict(120.8, {'GLY': 120.8}),
    ('CA', 'C', 'OXT'): DefaultValueDict(117.0),
}

BOND_ANGLE_SIDECHAIN = {
    'CYS': {('CA', 'CB', 'SG'): 114.4},
    'ASP': {('CA', 'CB', 'CG'): 112.6, ('CB', 'CG', 'OD1'): 118.4, 
        ('CB', 'CG', 'OD2'): 118.4},   
    'GLU': {('CA', 'CB', 'CG'): 114.1, ('CB', 'CG', 'CD'): 112.6,
        ('CG', 'CD', 'OE1'): 118.4, ('CG', 'CD', 'OE2'): 118.4},
    'PHE': {('CA', 'CB', 'CG'): 113.8, ('CB', 'CG', 'CD1'): 120.7, 
        ('CD1', 'CG', 'CD2'): 118.6, ('CG', 'CD1', 'CE1'): 120.7,
        ('CG', 'CD2', 'CE2'): 120.7, ('CD1', 'CE1', 'CZ'): 120.0,
        ('CD2', 'CE2', 'CZ'): 120.0, ('CE2', 'CZ', 'CE1'): 120.0},
    'HIS': {('CA', 'CB', 'CG'): 113.8, ('CB', 'CG', 'ND1'): 121.6,
        ('ND1', 'CG', 'CD2'): 109.3, ('CG', 'ND1', 'CE1'): 105.6,
        ('CG', 'CD2', 'NE2'): 106.5, ('ND1', 'CE1', 'NE2'): 111.7,
        ('CD2', 'NE2', 'CE1'): 107.0},
    'ILE': {('CA', 'CB', 'CG1'): 110.4, ('CA', 'CB', 'CG2'): 110.5, 
        ('CB', 'CG1', 'CD1'): 113.8}, 
    'LYS': {('CA', 'CB', 'CG'): 114.1, ('CB', 'CG', 'CD'): 111.3,
        ('CG', 'CD', 'CE'): 111.3, ('CD', 'CE', 'NZ'): 111.9},
    'LEU': {('CA', 'CB', 'CG'): 116.3, ('CB', 'CG', 'CD1'): 110.7, 
        ('CB', 'CG', 'CD2'): 110.7}, 
    'MET': {('CA', 'CB', 'CG'): 114.1, ('CB', 'CG', 'SD'): 112.7,
        ('CG', 'SD', 'CE'): 100.9},
    'ASN': {('CA', 'CB', 'CG'): 112.6, ('CB', 'CG', 'OD1'): 120.8,
        ('CB', 'CG', 'ND2'): 116.4},
    'PRO': {('CA', 'CB', 'CG'): 104.5, ('CB', 'CG', 'CD'): 106.1},
    'GLN': {('CA', 'CB', 'CG'): 114.1, ('CB', 'CG', 'CD'): 112.6,
        ('CG', 'CD', 'OE1'): 120.8, ('CG', 'CD', 'NE2'): 116.4},
    'ARG': {('CA', 'CB', 'CG'): 114.1, ('CB', 'CG', 'CD'): 111.3,
        ('CG', 'CD', 'NE'): 112.0, ('CD', 'NE', 'CZ'): 124.2,
        ('NE', 'CZ', 'NH1'): 120.0, ('NE', 'CZ', 'NH2'): 120.0},
    'SER': {('CA', 'CB', 'OG'): 111.1},
    'THR': {('CA', 'CB', 'OG1'): 109.6, ('CA', 'CB', 'CG2'): 110.5},
    'VAL': {('CA', 'CB', 'CG1'): 110.5, ('CA', 'CB', 'CG2'): 110.5},
    'TRP': {('CA', 'CB', 'CG'): 113.6, ('CB', 'CG', 'CD1'): 126.8,
        ('CD1', 'CG', 'CD2'): 106.3, ('CG', 'CD1', 'NE1'): 110.2,
        ('CG', 'CD2', 'CE2'): 107.2, ('CE2', 'CD2', 'CE3'): 118.9,
        ('CD2', 'CE2', 'CZ2'): 122.4, ('CD2', 'CE', 'CZ3'): 118.6,
        ('CE3', 'CZ3', 'CH2'): 121.1},
    'TYR': {( 'CA', 'CB', 'CG'): 113.9, ('CB', 'CG', 'CD1'): 120.8, 
        ('CD1', 'CG', 'CD2'): 118.4, ('CG', 'CD1', 'CE1'): 121.2,
        ('CG', 'CD2', 'CE2'): 121.2, ('CD1', 'CE1', 'CZ'): 119.6,
        ('CD2', 'CE2', 'CZ'): 119.6, ('CE1', 'CZ', 'OH'): 119.9},
}

BOND_ANGLE_PSEUDO_SIDECHAIN = {'ALA': 89.8283, 'CYS': 118.8, 'ASP': 117.5, 
    'GLU': 108.9, 'PHE': 126.1, 'GLY': 89.8283, 'HIS': 126.1, 'ILE': 100.3,
    'LYS': 111.7, 'LEU': 108.9, 'MET': 103.1, 'ASN': 118.8, 'PRO': 126.1,
    'GLN': 106.0, 'ARG': 114.6, 'SER': 91.7, 'THR': 94.5, 'VAL': 88.8,
    'TRP': 108.9, 'TYR': 126.1}



def get_bond_length(atom1, atom2, residue):
    """ Get the length of a the bond between two given atoms, according to the 
    residue.
    
    @param atom1: The first atom in the bond.
    @type atom1: str (one of the atoms defined in the L{ATOMS} constant)
    
    @param atom2: The second atom in the bond.
    @type atom2: str (one of the atoms defined in the L{ATOMS} constant)
    
    @rtype: float
    @return: The bond length. 
        
    """
    bond_length =  0.0    

    # Hydrogen atoms are not from Engh, Huber, 1991.
    if ATOMS.index(atom2) >= ATOMS.index('H') and \
        ATOMS.index(atom2) <= ATOMS.index('H3'):
        bond_length = 1.0
    else:
        # atoms are ordered in the bond length dictionary to avoid
        # duplication of entries
        if atom1 > atom2:
            atom1, atom2 = atom2, atom1        
        
        try:
            bond_length = BOND_LENGTH_MAINCHAIN[(atom1, atom2)][residue]
        except KeyError:
            try:
                bond_length = BOND_LENGTH_SIDECHAIN[residue][atom1, atom2]
            except KeyError:
                    bond_length = BOND_LENGTH_PSEUDO_SIDECHAIN[residue]       
            
    assert bond_length != 0, "Bond length could not be obtained."
    return bond_length
        
        
def get_bond_angle(atom1, atom2, atom3, residue):
    """ Get the angle of the bond between three atoms, according to the residue.
    
    @param atom1: The first atom in the bond.
    @type atom1: str (one of the atoms defined in the L{ATOMS} constant)
    
    @param atom2: The second atom in the bond.
    @type atom2: str (one of the atoms defined in the L{ATOMS} constant)
    
    @param atom3: The third atom in the bond.
    @type atom3: str (one of the atoms defined in the L{ATOMS} constant)
    
    @rtype: float
    @return: The bond angle. 
    
    """
    bond_angle = 0.0

    try:
        bond_angle = BOND_ANGLE_MAINCHAIN[(atom1, atom2, atom3)][residue]
    except KeyError:
        try:
            bond_angle = BOND_ANGLE_SIDECHAIN[residue][atom1, atom2, atom3]
        except KeyError:
                bond_angle = BOND_ANGLE_PSEUDO_SIDECHAIN[residue]       
            
    assert bond_angle != 0, "Bond angle could not be obtained."
        
    return bond_angle * math.pi / 180.0
    


