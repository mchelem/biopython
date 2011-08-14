try:
    import mocapy
except ImportError:        
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Mocapy was not found. "
        "Install the Mocapy library if you "
        "want to use Bio.PDB.RNAStructure.")

from Bio.PDB.TorusDBN.TorusDBNModel import TorusDBNModel
from Bio.PDB.TorusDBN.TorusDBNTrainer import TorusDBNTrainer

