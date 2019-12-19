"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

The script contains definitions of custom Exceptions
for GRAFIMO

"""

import sys
from shutil import which
import numpy as np

"""
    definition of constant variables
"""

DNA_ALPHABET=['A', 'C', 'G', 'T'] # dna alphabet (we ignore the N and the IUPAC symbols)
REV_COMPL={'A':'T', 'C':'G', 'G':'C', 'T':'A'}
LOG_FACTOR=1.44269504
RANGE=1000
CHROMS_LIST=[str(i) for i in range(1, 23)] + ['X', 'Y']
EXT_DEPS = ['tabix', 'vg']

"""
    functions from utils.py
"""
def die(code):
    """
        Exit from the execution
        ----
        Parameters:
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
        Parameters:
            lst1 (list) : first list
            lst2 (list) : second list
        ----
        Returns:
            (bool)
    """

    if len(lst1)==len(lst2) and set(lst1)==set(lst2):
        return True

    return False


def check_deps():
    """
        Checks if all external dependencies needed by
        GRAFIMO (vg and tabix) are satisfied
        ----
        Parameters:
            None
        ----
        Returns:
            sat (bool) : set to False if at least one 
                            dependency is not satisfied
            deps_not_sats (list) : list containing the 
                                    dependencies that are
                                    not satisfied
    """

    deps_not_sats = []
    sat = True

    for dep in EXT_DEPS:
        if not which(dep) is not None:
            deps_not_sats.append(dep)
            sat = False

    return sat, deps_not_sats


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
        Parameters:
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
        Parameters:
            value (np.double) : value of which the log2 will be computed
        ----
        Returns:
            (np.double)
    """

    return (np.log(value)*LOG_FACTOR)


def correct_path(path, path_id='', file_format=''):

    if path[-1:] == '/':
        new_path = ''.join([path, path_id, file_format])
    else:
        new_path = ''.join([path, '/', path_id, file_format])

    return new_path


def unique_lst(lst):
    """
        Get the unique values inside a list
        ----
        Parameters:
            lst (list) : list of values
        ----
        Returns:
            unique_lst (list) : list of the unique values in lst
    """

    assert(len(lst) > 0)

    unique_lst = []

    for el in lst:

        if el not in unique_lst:
            unique_lst.append(el)

    assert(len(unique_lst) > 0)

    return unique_lst


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
        Print the progress bar in the sequence scoring process and graph extraction
        process
        ----
        Parameters:
            iteration (int)
            total (int)
            prefix (str)
            suffix (str)
            decimals (int)
            length (int)
            fill (str)
            printEnd (str)
        ----
        Returns
            None
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()   
