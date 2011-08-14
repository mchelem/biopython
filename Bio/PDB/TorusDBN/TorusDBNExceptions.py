# Copyright (C) 2011 by Michele Silva (michele.silva@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from Bio import BiopythonWarning


class TorusDBNException(Exception):
    """
    General TorusDBN error.
    """
    
    
class TorusDBNBuildPolypeptideException(TorusDBNException):
    """
    Can't build a list of Polypeptide objects.
    """
