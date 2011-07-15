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

from Bio.PDB import Vector
from numpy import array, matrix
from Bio.PDB.Vector import calc_angle
from math import sqrt, cos, sin, pi

class ConversionException(Exception):
    def __init__(self, value, position=None):
        self.value = value
        self.position = position
        
    def __str__(self):
        return ("ReconstructException:" + repr(self.value))

    def get_position(self):
        return self.position


def fromDihedral(A, B, C, bond_BC, bond_CD, angle_BCD, torsion_BC):
    """
    Calculate a new point from three points A,B,C, an dihedral angle
    (torsion_BC) and bond angle (angle_BCD). Done by putting a
    coordinate system in C and then going from spherical coordinates
    and then translating the coordinate system to C.

    See Parsons et al. for details
    """

    # Turn the bond angle into an angle in [pi/2, pi]
    angle_BCD = pi - angle_BCD 

    # Calculate position of D from spherical coordinate representation
    D2 = Vector(cos(angle_BCD), cos(torsion_BC)*sin(angle_BCD), sin(torsion_BC)*sin(angle_BCD))**bond_CD

    # Calculate rotation matrix M
    bc = (C-B) / float(bond_BC)  # Normalized by previous bond length
    n = (B-A)**bc
    n.normalize()                # Normalized by calculation
    nXbc = n**bc                 # Normalized by definition

    M = (array([bc.get_array(), nXbc.get_array(), n.get_array()])).T

    # Calculate position of D by rotation and translation
    D = D2.left_multiply(M) + C
    
    return D

def from3atoms(A,B,C,alpha,beta,l,choice):
    gamma = calc_angle(B,A,C)

    Dx = l*cos(alpha)
    Dy = l*(cos(beta)-cos(alpha)*cos(gamma))/sin(gamma)
    Dz = sqrt(l**2-Dx**2-Dy**2)

    Dab   = Dx-cos(gamma)/sin(gamma)*Dy
    Dac   = Dy/sin(gamma)
    Dabac = Dz

    AB    = (B-A).normalized()
    AC    = (C-A).normalized()
    ABxAC = (AB**AC).normalized()

    if choice == 0:
        AD = AB**Dab + AC**Dac + ABxAC**Dabac
    else:
        AD = AB**Dab + AC**Dac - ABxAC**Dabac

    return A+AD


def from4atoms(OA, OB, OC, OD, u, v, l, m, accept_no_solution=True):
    u = pi - u
    v = pi - v
    
    a = OA - OC
    b = OB - OD
    
    ddiff = OC - OD + a - b # (BA) = (B-A)
    x = a.norm() * l * cos(u)
    y = b.norm() * m * cos(v)

    rx = x
    ry = y - b * ddiff

    r1 = a
    r2 = b

    if(abs(r1[0]) < 0.000001):
        if(abs(r2[0]) < 0.000001):
            print "WARNING: r1 = r2 = 0.0"        
        # Swap rows
        tmp = r1
        r1 = r2
        r2 = tmp
        tmp = rx
        ry = rx
        rx = tmp

    # Reduce the matrix
    factor = r1[0]
    r1 = r1 / factor
    rx = rx / factor    
    factor = r2[0]
    r2 = r2 - (r1 ** factor)
    ry = ry - rx * factor
    factor = r2[1]
    r2 = r2 / factor
    ry = ry / factor
    factor = r1[1]
    r1 = r1 - r2 ** factor
    rx = rx - ry * factor

    # Make solution space vectors
    alpha = r1[2]
    beta  = r2[2]
    gamma = rx
    delta = ry
    u = Vector(-alpha,-beta,1.0)
    v = Vector(gamma, delta, 0)

    # Solve quadratic equation for norm of c
    acoef = (u.norm())**2
    bcoef = 2*(u*v)
    ccoef = (v.norm())**2 - l**2
    disc = bcoef**2 - 4*acoef*ccoef

    if(disc < 0.0):
        # This is the sick case where we can't find _any_ solution
        if not accept_no_solution:
            raise ValueError, "from4atoms: no solution found (disc=%f)" % disc
        else:
            disc = 0
    
    x1 = (-bcoef - sqrt(disc))/(2*acoef)
    x2 = (-bcoef + sqrt(disc))/(2*acoef)

    # Create the two c-vectors
    c1 = u ** x1 + v
    c2 = u ** x2 + v

    # The two candidate E-poitns
    E1 = OA + c1
    E2 = OA + c2

    # The two candidate d-vectors
    d1 = E1 - OB
    d2 = E2 - OB

    # Pick the one with smallest norm difference
    d1norm = d1.norm()
    d2norm = d2.norm()
    diff1 = abs(d1norm - m)
    diff2 = abs(d2norm - m)

    if(diff1 < diff2):
        return E1
    else:
        return E2


def place_base(N_B, K_B, C_B, N_Bp, K_Bp, C_Bp, transform_list):
    # Calculate the basis for B'' and basis change matrix from B' to B''
    Xpp_Bp = (K_Bp-N_Bp).normalized()
    Zpp_Bp = ((K_Bp-N_Bp)**(C_Bp-N_Bp)).normalized()
    Ypp_Bp = Zpp_Bp**Xpp_Bp

    Bp_to_Bpp = ((matrix([Xpp_Bp.get_array(), Ypp_Bp.get_array(), Zpp_Bp.get_array()])).T).I

    # Calculate the basis vectors of B'' in B and basis change matrix from B'' to B
    Xpp_B = (K_B-N_B).normalized()
    Zpp_B = ((K_B-N_B)**(C_B-N_B)).normalized()
    Ypp_B = Zpp_B**Xpp_B

    Bpp_to_B = (matrix([Xpp_B.get_array(), Ypp_B.get_array(), Zpp_B.get_array()])).T

    # Calculate the change matrix from B' to B
    Bp_to_B = (Bpp_to_B * Bp_to_Bpp).A
    
    # Do the transformation for the whole list
    transform_list = map(lambda Vp: (Vp-N_Bp).left_multiply(Bp_to_B)+N_B, transform_list)
    
    return transform_list
