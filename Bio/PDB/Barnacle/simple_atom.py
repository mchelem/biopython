# Copyright (C) 2008 by Jes Frellsen, Ida Moltke and Martin Thiim
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""
This module contains a function that simplifies creating an atom in Bio.PDB.

From Biopython version 1.53 the PDB column element was introduced. To
accommodate this, two version of the sAtom function has be
written. One from pre 1.53 Bio.PDB and one for post 1.53.

The version of Atom class is decided by checking for the 'element'
argument in the Atom constructor.

"""

from Bio.PDB import Atom
import inspect


def full_atom_name(name):
    if len(name) <= 3:
        return " %-3s" % name
    elif len(name) == 4:
        return name
    else:
        raise Exception, "Atom names can be no longer than 4 characters!"


def sAtom_pre_1_53(name, vector):
    """
    A wrapper function for creating a simple atom for Bio.PDB pre version 1.53.
    """
    return Atom.Atom(name, vector.get_array(), 0, 1, " ", full_atom_name(name), None)


def sAtom_post_1_53(name, vector):
    """
    A wrapper function for creating a simple atom for Bio.PDB post version 1.53.
    """
    full_name = full_atom_name(name)
    atom = Atom.Atom(name, vector.get_array(), 0, 1, " ", full_name, None, name[0])
    return atom

# Define the sAtom function depending on the Bio.PDB version
if 'element' in inspect.getargspec(Atom.Atom.__init__)[0]:
    sAtom = sAtom_post_1_53
else:
    sAtom = sAtom_pre_1_53

