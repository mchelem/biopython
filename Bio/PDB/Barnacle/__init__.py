# Copyright (C) 2008 by Jes Frellsen, Ida Moltke and Martin Thiim
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

try:
    import mocapy
except ImportError:        
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Mocapy was not found. "
        "Install the Mocapy library if you "
        "want to use Bio.PDB.Barnacle.")
