"""Functions and constant variables used in GRAFIMO code"""


from grafimo.GRAFIMOException import NoDataFrameException
from typing import List, Optional, Tuple
from shutil import which
import pandas as pd
import numpy as np
import sys
import os


#-----------------------------------------------------------------------
# constant vars
#-----------------------------------------------------------------------
DNA_ALPHABET = ['A', 'C', 'G', 'T'] 
REV_COMPL = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
PSEUDObg = np.double(0.0000005)
LOG_FACTOR = 1.44269504
RANGE = 1000
CHROMS_LIST = [str(i) for i in range(1, 23)] + ['X', 'Y']
DEFAULT_OUTDIR = "default_out_dir_name"
EXT_DEPS = ['tabix', 'vg', 'dot']
SOURCE = 'grafimo'
TP = 'nucleotide_motif'
PHASE = '.'


#-----------------------------------------------------------------------
# functions
#-----------------------------------------------------------------------
def die(code: int) -> None:
    """Stop the execution and exit.

    Parameters
    ----------
    code : int
        stop code
    """

    sys.exit(code)

# end of die()


def sigint_handler() -> None:
    """Print a message when a SIGINT is caught and exit."""

    print("\nCaught SIGINT. GRAFIMO will exit")
    die(2)

# end of sigint_handler()


def isListEqual(lst1: List, lst2: List) -> bool:
    """Check if two lists contain the same distinct elements.

    Parameters
    ----------
    lst1 : list
        list 1
    lst2 : list
        list 2
    
    Returns
    -------
    bool
    """

    if (len(lst1) == len(lst2) and set(lst1) == set(lst2)):
        return True

    return False

# end of isListEqual()


def initialize_chroms_list(args_chroms: List[str]):
    """Initialize a list of chromomsome to use during GRAFIMO analysis.

    Parameters
    ----------
    args_chroms : list
        list of user chromosomes
    
    Returns
    -------
    list
        chromosome list
    """

    if not args_chroms:
        chroms = CHROMS_LIST # all the chromosomes

    else:
        chroms = args_chroms # the given chromosomes

    return chroms

# end of initialize_chroms_list()


def check_deps() -> Tuple[bool, List[str]]:
    """Check if the external dependencies needed to run GRAFIMO are 
    satisfied.

    Returns
    -------
    bool
        check result
    list
        list of dependencies not satisfied
    """

    deps_not_sats: List[str] = list()
    sat: bool = True

    for dep in EXT_DEPS:
        if not which(dep) is not None:
            deps_not_sats.append(dep)
            sat = False

    return sat, deps_not_sats

# end of check_deps()


def isJaspar_ff(motif_file: str) -> bool:
    """Check if the given motif file is a JASPAR file
        
    Parameters
    ----------
    motif_file : str 
        motif file
        
    Returns
    -------
    bool
        check result
    """

    if not os.path.isfile(motif_file):
        errmsg: str = "\n\nERROR: unable to find %s" % motif_file
        raise FileNotFoundError(errmsg)

    if motif_file and isinstance(motif_file, str):
        ff: str = motif_file.split('.')[-1]

        if ff == 'jaspar':
            return True
        else:
            return False

    else:
        return False

# end of isJaspar_ff()


def isMEME_ff(motif_file: str) -> bool:
    """Check if the given motif file is a MEME file
        
    Parameters
    ----------
    motif_file : str 
        motif file
        
    Returns
    -------
    bool
        check result
    """

    if not os.path.isfile(motif_file):
        errmsg: str = "\n\nERROR: unable to find %s" % motif_file
        raise FileNotFoundError(errmsg)

    if motif_file and isinstance(motif_file, str):
        ff: str = motif_file.split('.')[-1]

        if ff == 'meme':
            return True
        else:
            return False

    else:
        return False 

# end of isMEME_ff()


def almost_equal(value1: np.double,
                 value2: np.double,
                 slope: np.double
) -> bool:
    """Check if two values are 'close' to each other. given a degree
    of tolerance
    
    Parameters
    ----------
    value1 : np.double
        first value
    value2 : np.double 
        second value
    slope : np.double
        tolerance
        
    Returns
    -------
    bool
    """

    if ((value1 - slope) > value2 or (value1 + slope) < value2):
        return False
    
    return True

# end of almost_equal()


def lg2(value: np.double) -> np.double:
    """C like computation of log2
        
    Parameters
    ----------
    value : np.double 
        value
        
    Returns
    np.double
        log2 of the input value
    """

    log2value = (np.log(value) * LOG_FACTOR)
    return log2value

# end of lg2()


def unique_lst(lst: List, size: Optional[int] = None) -> List:
    """Get the unique values contained in a list and store them in 
    another one.

    Parameters
    ----------
    lst : list
        list of values
    size : int, optional
        number of unique values that the list has to contain

    Returns
    -------
    list
        list of the unique values
    """

    assert (len(lst) > 0)

    unique_lst: List = list() 
    el_num: int = 0

    for el in lst:

        if el not in unique_lst:
            unique_lst.append(el)
            el_num += 1

        if size != None and el_num == size:  # size limit reached
            break

    assert (len(unique_lst) > 0)

    return unique_lst

# end of unique_lst()


def list_data(data: pd.DataFrame, qvalue: bool) -> List:
    """Convert a pandas DataFrame in a list containign each dataframe 
    column as a list of values

    Parameters
    ----------
    data : pandas.DataFrame
        input DataFrame
    qvalue : bool
        if True the column of q-values has to be considered

    Returns
    -------
    list
        list containing DataFrame's columns as list of values
    """


    if not isinstance(data, pd.DataFrame):
        errmsg: str = "\n\nERROR: not allowed data type given"
        raise NoDataFrameException(errmsg)

    assert len(data.columns) <= 12
    assert len(data.columns) >= 11

    seqnames: List[str] = data['sequence_name'].to_list()
    starts: List[int] = data['start'].to_list()
    stops: List[int] = data['stop'].to_list()
    scores: List[np.double] = data['score'].to_list()
    strands: List[str] = data['strand'].to_list()
    motifIDs: List[str] = data['motif_id'].to_list()
    motifNames: List[str] = data['motif_alt_id'].to_list()
    pvalues: List[np.double] = data['p-value'].to_list()
    sequences: List[str] = data['matched_sequence'].to_list()
    frequencies:List[int] = data['haplotype_frequency'].to_list()
    references: List[str] = data['reference'].to_list()

    if qvalue:
        qvalues: List[np.double] = data['q-value'].to_list()

    if qvalue:
        summary = [motifIDs, motifNames, seqnames, starts, stops, strands, scores,
                   pvalues, sequences, frequencies, references, qvalues]
    else:
        summary = [motifIDs, motifNames, seqnames, starts, stops, strands, scores,
                   pvalues, sequences, frequencies, references]

    summary_len: int = len(motifIDs)

    assert summary_len == len(data.index)
    assert summary_len == len(motifNames)
    assert summary_len == len(seqnames)
    assert summary_len == len(starts)
    assert summary_len == len(stops)
    assert summary_len == len(strands)
    assert summary_len == len(scores)
    assert summary_len == len(pvalues)
    assert summary_len == len(sequences)
    assert summary_len == len(frequencies)
    assert summary_len == len(references)

    if qvalue:
        assert summary_len == len(qvalues)

    return summary

# end of list_data()


def printProgressBar(iteration: int,
                     total: int,
                     prefix: Optional[str] = '',
                     suffix: Optional[str] = '',
                     decimals: Optional[int] = 1,
                     length: Optional[int] = 50,
                     fill: Optional[str] = '=',
                     printEnd: Optional[str] ="\r"
) -> None:
    """Print a progress bar.

    The progress bar is printed while processing MEME files conatining 
    many motifs, while building VGs, while extracting motif candidate 
    and while scoring the motif candidates.

    Parameters
    ----------
    iteration : int
        iteration number
    total : int
        total number of iterations to do
    prefix : str
        string to print in front of the bar
    suffix : str
        string t print at the end of the bar
    decimals : int
        number of decmal digits to display
    length : int
        length of the bar (# characters to use)
    fill : str
        string to fill the bar
    printEnd : str
        string to print at end of the whole bar 'structure'
    """

    percent: str = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength: int = int(length * iteration // total)
    # "allocate space" for the bar
    bar: float = fill * filledLength + ' ' * (length - filledLength)  

    print('\r%s [%s] %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)

    # new line when the bar is completely filled
    if iteration == total:
        print()

# end of printProgressBar()

