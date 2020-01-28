"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Utilities used in the different steps of GRAFIMO
pipeline

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
EXT_DEPS = ['tabix', 'vg', 'dot']


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


def initialize_chroms_list(args_chroms, chroms_lst = CHROMS_LIST):
    """
        Initialize the list of chromosomes that will be 
        considered by GRAFIMO, during its run.
        ----
        Parameters:
            args_chroms (list) : chromosome list as it is obtained 
                                    from the user input
            chroms_lst (list) : list of all the chromosome 
        ----
        Returns:

    """

    chroms = [] # result variable

    if not args_chroms:
        chroms = chroms_lst
    
    else:
        chroms = args_chroms

    return chroms


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


def unique_lst(lst, size = None):
    """
        Get the unique values inside a list
        ----
        Parameters:
            lst (list) : list of values
            size (int) : number of unique elements the list must contain.
                         By default it is set to None, then all the element
                         of the list will be checked.
        ----
        Returns:
            unique_lst (list) : list of the unique values in lst
    """

    assert(len(lst) > 0)

    unique_lst = []
    el_num = 0

    for el in lst:

        if el not in unique_lst:
            unique_lst.append(el)
            el_num += 1

        if size != None and el_num == size: # size limit reached
            break

    assert(len(unique_lst) > 0)

    return unique_lst


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 50, fill = '=', printEnd = "\r"):
    """
        Print the progress bar in the sequence scoring process and graph extraction
        process
        ----
        Parameters:
            iteration (int) : fraction of work done
            total (int) : total amount of work to do
            prefix (str) : string to put in fornt of the bar
            suffix (str) : string to put at the end of the bar
            decimals (int) : number of digits after the dot, in the
                                representation of the percentage of work done
            length (int) : length of the bar
            fill (str) : character with will be filled the bar
            printEnd (str) : what will be done when the bar has been 
                                completely printed
        ----
        Returns:
            None
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + ' ' * (length - filledLength) # "allocate space" for the bar
    
    # print the bar
    print('\r%s [%s] %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    
    # new line when the bar is completely filled
    if iteration == total:
        print()   

