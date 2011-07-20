try:
    import mocapy
except ImportError:        
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Mocapy was not found. "
        "Install the Mocapy library if you "
        "want to use Bio.PDB.RNAStructure.")
