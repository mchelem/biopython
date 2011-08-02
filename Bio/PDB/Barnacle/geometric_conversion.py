# Copyright (C) 2008 by Jes Frellsen, Ida Moltke and Martin Thiim
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from Bio.PDB import Vector, PDBIO
from Bio.PDB.Residue import Residue
from Bio.PDB.StructureBuilder import StructureBuilder
from math import cos, sin, pi
from NucleicAcid import has_base_type, is_purine, get_one_letter_code, get_three_letter_code
from Bio.PDB.Barnacle.geometric_constants import *
from Bio.PDB.Barnacle.geometric_conversion_linalg import fromDihedral, from3atoms, from4atoms, place_base, ConversionException
from Bio.PDB.Barnacle.simple_atom import sAtom

def atoms_to_residues(atoms, atomsPerRes, seqid, resnames=None):
    """
    Builds a list of residues from a list of atoms, where each residue
    i the list has 'atomsPerRes' atoms per residue. All residues are
    set to 'A', if no residue names are given.
    """

    # Check the arguments
    if (len(atoms) % atomsPerRes) != 0:
        raise ValueError, "Mismatch between list length and atoms per residues."

    if resnames==None:
        resnames = ["  A" for _ in xrange(len(atoms)/atomsPerRes)]


    residues = []
    res = None
    
    for (i,atom) in enumerate(atoms):
        # Residue number
        resnumber = i // atomsPerRes

        # Create a new residue
        if (i % atomsPerRes) == 0:
            if res!=None:
                residues.append(res)
            res = Residue((" ", resnumber, " "), resnames[resnumber], seqid)
            res.add(atom)
        # Add an atom
        else:
            res.add(atom)

    if res!=None:
        residues.append(res)

    return residues


def residues_to_struct(name_residues_list, structID):
    """
    Build a structure from a list of (chain name, residues) - each as
    a chain.
    """
    
    builder = StructureBuilder()
    builder.init_structure(structID)
    builder.init_model(0)

    for (name, residues) in name_residues_list:
        builder.init_chain(name)
        for res in residues:
            builder.chain.add(res)

    return builder.get_structure()


def convert_from_bb_dihedrals(dihedrals_list, check_valid_dihedrals=True):
    """
    Construct the RNA backbone from a list of tuples (alpha, beta,
    gamma, delta, epsilon, zeta) of backbone dihedrals
    """ 
    if check_valid_dihedrals:
        # Check that the list is correct (only the first tuple may start with none and the last may end with two nones)
        if None in dihedrals_list[0][1:] or None in dihedrals_list[-1][:-2] or \
               (len(dihedrals_list)>2 and filter(lambda tup: None in tup, dihedrals_list[1:-1]) != []):
            raise ValueError, "The list of dihedrals is not valid"

    # List of Residues that is the result
    atoms = []

    # Setup the angles and distances to use
    atomnames = bb_atomnames

    # Use different angles and distances based on delta [3]
    if dihedrals_list[0][3] < 1.8:
        distances = bb_distances_3e
        angles = bb_angles_3e
    else:
        distances = bb_distances_2e
        angles = bb_angles_2e

    # Start in the origin and place the two next atoms in the x-y-plane
    P  = sAtom(atomnames[0], Vector(0,0,0))
    O5 = sAtom(atomnames[1], Vector(distances[0], 0, 0))
    C5 = sAtom(atomnames[2], Vector(distances[0], 0, 0) + Vector(cos(pi-angles[1]), sin(pi-angles[1]), 0)**distances[1])

    atoms = [P, O5, C5]

    n = len(atomnames)

    # Loop invariant:
    #  - Start: the list consists of whole sets of P-O5-C
    #  - During: last pair is CA-C-N
    #  - End: the list consists of whole sets of N-CA-C
    for (i, dihedrals) in enumerate(dihedrals_list):
        
        # Make sure that dihedrals[0][0], dihedrals[-1][-2] and dihedrals[-1][-1] are None
        if i==0 and dihedrals[0]!=None:
            dihedrals = list(dihedrals)
            dihedrals[0] = None
        if i==len(dihedrals_list)-1 and (dihedrals[-2]!=None or dihedrals[-1]!=None):
            dihedrals = list(dihedrals)
            dihedrals[-2] = dihedrals[-1] = None

        # Use different angles and distances based on delta [3]
        if dihedrals[3] < 1.8:
            distances = bb_distances_3e
            angles = bb_angles_3e
        else:
            distances = bb_distances_2e
            angles = bb_angles_2e

        # Place the atoms
        for (j, dihedral) in enumerate(dihedrals):
            if dihedral!=None:
                A, B, C = atoms[-3:]
                D = fromDihedral(A.get_vector(), B.get_vector(), C.get_vector(), \
                                 distances[j%n], distances[(j+1)%n], angles[(j+1)%n], dihedral)
                atoms.append(sAtom(atomnames[(j+2)%n], D))

    return atoms


def convert_from_dihedrals(dihedrals_list, sequence=None, bb_atoms=None, accept_no_solution=True,
                           add_sugar=True, add_base=True):
    """
    Construct a list of residues from a list of tuples of dihedral angles in the order:
    [(alpha, beta, gamma, delta, epsilon, zeta, chi), (alpha, beta ..., zeta, chi), ...]

    The angles must be in radians.

    @param dihedral_list: the list of tuples of dihedral angles.
    @type dihedral_list: L{(float, float, float, float, float, float, float)}

    @param sequence: optional string of nucleic acids. If set to None
                     all residues will be of type A
    @type sequence: string or None

    @param bb_atoms: optional list of back bone atoms. If this list is
                     given, the back bone will not be constructed from
                     the dihedral angles given.
    @type bb_atoms: L{Atom} or None
    
    @param accept_no_solution: If this option is set to false, the
                               function will raise an exception if the
                               method fails to find a solution for the
                               coordinates of C1'. This could indicate
                               that the dihedral angles given as
                               argument are 'bad'. If the options is
                               set to true, no exception will be
                               raised.
    @type accept_no_solution: bool

    @param add_sugar: indicated if the atoms in the sugar ring should be placed.
    @type add_sugar: bool

    @param add_base: indicated if the atoms in the base should be
                     placed. Ineffective if add_sugar is false
    @type add_base: bool

    @return: A list of residues
    @rtype: L{Resudue}
    
    """ 
    
    # Construct the backbone
    if bb_atoms==None:
        bb_dihedrals = map(lambda tup: tup[:6], dihedrals_list)
        bb_atoms = convert_from_bb_dihedrals(bb_dihedrals)

    if sequence!=None:
        sequence = map(get_three_letter_code, list(sequence))

    residues = atoms_to_residues(bb_atoms, 6, "0", resnames=sequence)
   
    # Construct each sugar and base
    for (resnumber, res) in enumerate(residues):

        # Attached oxygens OP1 and OP2 to the backbone P
        if resnumber!=0:
            O3p_prev = residues[resnumber-1]["O3'"].get_vector()
            P   = res["P"].get_vector()
            O5p = res["O5'"].get_vector()

            O1P = from3atoms(P, O3p_prev, O5p, a_O3pPO1P, a_O5pPO1P, d_PO1P, 0)
            O2P = from3atoms(P, O3p_prev, O5p, a_O3pPO2P, a_O5pPO2P, d_PO2P, 1)

            res.add(sAtom("OP1", O1P))
            res.add(sAtom("OP2", O2P))

        # Place the atoms in the sugar ring
        if add_sugar:
            
            # Make reference to the vectors of atom that will be used
            C5p = res["C5'"].get_vector()
            C4p = res["C4'"].get_vector()
            C3p = res["C3'"].get_vector()
            O3p = res["O3'"].get_vector()

            # Use different angles and distances based on delta [3]
            if dihedrals_list[resnumber][3] < 1.8:
                a_C5pC4pO4p = a_3e_C5pC4pO4p
                a_C3pC4pO4p = a_3e_C3pC4pO4p
                a_C2pC3pC4p = a_3e_C2pC3pC4p
                a_C2pC3pO3p = a_3e_C2pC3pO3p
                a_C4pO4pC1p = a_3e_C4pO4pC1p
                a_C1pC2pC3p = a_3e_C1pC2pC3p
                a_O4pC1pN   = a_3e_O4pC1pN
                a_C2pC1pN   = a_3e_C2pC1pN
                d_C4pO4p    = d_3e_C4pO4p
                d_C2pC3p    = d_3e_C2pC3p
                d_O4pC1p    = d_3e_O4pC1p
                d_C1pC2p    = d_3e_C1pC2p
                d_C1pN      = d_3e_C1pN

                a_C1pC2pO2p = a_3e_C1pC2pO2p
                a_C3pC2pO2p = a_3e_C3pC2pO2p
                d_C2pO2p    = d_3e_C2pO2p


            else:
                a_C5pC4pO4p = a_2e_C5pC4pO4p
                a_C3pC4pO4p = a_2e_C3pC4pO4p
                a_C2pC3pC4p = a_2e_C2pC3pC4p
                a_C2pC3pO3p = a_2e_C2pC3pO3p
                a_C4pO4pC1p = a_2e_C4pO4pC1p
                a_C1pC2pC3p = a_2e_C1pC2pC3p
                a_O4pC1pN   = a_2e_O4pC1pN
                a_C2pC1pN   = a_2e_C2pC1pN
                d_C4pO4p    = d_2e_C4pO4p
                d_C2pC3p    = d_2e_C2pC3p
                d_O4pC1p    = d_2e_O4pC1p
                d_C1pC2p    = d_2e_C1pC2p
                d_C1pN      = d_2e_C1pN

                a_C1pC2pO2p = a_2e_C1pC2pO2p
                a_C3pC2pO2p = a_2e_C3pC2pO2p
                d_C2pO2p    = d_2e_C2pO2p

            # From C5p,C4p and C3p place 04p
            O4p = from3atoms(C4p,C5p,C3p, a_C5pC4pO4p, a_C3pC4pO4p, d_C4pO4p,0) # Always solution 0
            res.add(sAtom("O4'", O4p))

            # From C4p,C3p and O3p place C2p
            C2p = from3atoms(C3p,C4p,O3p, a_C2pC3pC4p, a_C2pC3pO3p, d_C2pC3p,1) # Always solution 1
            res.add(sAtom("C2'", C2p))

            # From O4p and C2p place C1
            try:
                C1p = from4atoms(O4p,C2p,C4p,C3p,a_C4pO4pC1p,a_C1pC2pC3p,d_O4pC1p,d_C1pC2p, accept_no_solution=accept_no_solution)
            except ValueError, exception:
                raise ConversionException(0, resnumber)
            res.add(sAtom("C1'", C1p))

            # Attached oxygen O2' to C2' in the sugar ring
            O2p = from3atoms(C2p, C1p, C3p, a_C1pC2pO2p, a_C3pC2pO2p, d_C2pO2p, 1)
            res.add(sAtom("O2'", O2p))

            # Place the atoms in the base 
            if add_base:

                # From O4p,C1p and C2p place N
                N = from3atoms(C1p,O4p,C2p, a_O4pC1pN, a_C2pC1pN, d_C1pN,0)         # Always solution 0
                if is_purine(res):
                    res.add(sAtom("N9", N))
                else:
                    res.add(sAtom("N1", N))

                # From C2p,C1p,N and chi place C2/C4 and base
                if has_base_type("A", res):
                    d_NC = d_A_N9C4
                    a_CNC1p = a_A_C4N9C1p
                    base_dict, base_coords, base_names = A_dict, A_coords, A_names
                elif has_base_type("C", res):
                    d_NC = d_C_N1C2
                    a_CNC1p = a_C_C2N1C1p
                    base_dict = C_dict
                    base_dict, base_coords, base_names = C_dict, C_coords, C_names
                elif has_base_type("G",res):
                    d_NC = d_G_N9C4
                    a_CNC1p = a_G_C4N9C1p
                    base_dict, base_coords, base_names = G_dict, G_coords, G_names
                elif has_base_type("U",res):
                    d_NC = d_U_N1C2
                    a_CNC1p = a_U_C2N1C1p
                    base_dict, base_coords, base_names = U_dict, U_coords, U_names
                else:
                    print "Warning: Unsupported base, %s, met" % get_one_letter_code(res)

                chi = dihedrals_list[resnumber][6]
                if chi==None:
                    raise ValueError, "The list of chis is not valid"

                C24 = fromDihedral(O4p, C1p, N, d_C1pN, d_NC, a_CNC1p, chi)

                # Do the actual base placing
                if is_purine(res):
                    base = place_base(N, C24, C1p, base_dict["N9"], base_dict["C4"], base_dict["C1'"], base_coords[2:])
                else:
                    base = place_base(N, C24, C1p, base_dict["N1"], base_dict["C2"], base_dict["C1'"], base_coords[2:])

                base_atom_names = base_names[2:]

                for (i,atom) in enumerate(base):
                    res.add(sAtom(base_atom_names[i], atom))
  
    return residues

