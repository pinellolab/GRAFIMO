"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

The script contains definitions of custom Exceptions
for GRAFIMO

"""

import sys
import numpy as np

"""
    definition of constant variables
"""

DNA_ALPHABET=['A', 'C', 'G', 'T'] # dna alphabet (we ignore the N and the IUPAC symbols)
REV_COMPL={'A':'T', 'C':'G', 'G':'C', 'T':'A'}
LOG_FACTOR=1.44269504
RANGE=1000

"""
    functions from utils.py
"""
def die(code):
    """
        Exit from the execution
        ----
        Params:
            code (int) : exit code
        ----
        Returns:
             None
    """

    sys.exit(code)

def isListEqual(lst1, lst2):
    """
        Compare two lists if they are equal
        ----
        Params:
            lst1 (list) : first list
            lst2 (list) : second list
        ----
        Returns:
            (bool)
    """

    if len(lst1)==len(lst2) and set(lst1)==set(lst2):
        return True

    return False


def isJaspar_ff(motif_file):
    """
        Check if a file is in .jaspar format
        ----
        Parameters:
            motif_file (str) : path to the motif file
        ----
        Returns:
            (bool)
    """

    if motif_file and isinstance(motif_file, str):
        ff = motif_file.split('.')[-1]

        if ff == 'jaspar':
            return True
        else:
            return False

    else:
        return False  # the motif file was not given as a path or the path is of length 0


def isMEME_ff(motif_file):
    """
        Check if the given file is in .meme format
        ----
        Parameters:
            motif_file (str) : path to the motif file
        ----
        Returns:
            (bool)
    """

    if motif_file and isinstance(motif_file, str):
        ff = motif_file.split('.')[-1]

        if ff == 'meme':
            return True
        else:
            return False

    else:
        return False  # the motif file was not given or the path is empty


def almost_equal(value1, value2, slope):
    """
        Computes if two values are close considering a slope as degree
        of freedom
        ----
        Params:
            value1 (np.double) : first value
            value2 (np.double) : second value
            slope (np.double) : tolerance
        ----
        Returns:
             (bool)
    """

    if (value1 - slope) > value2 or (value1 + slope) < value2:
        return False
    else:
        return True

def lg2(value):
    """
        C-like implementation of the log2 with a faster running time
        ----
        Params:
            value (np.double) : value of which the log2 will be computed
        ----
        Returns:
            (np.double)
    """

    return (np.log(value)*LOG_FACTOR)

