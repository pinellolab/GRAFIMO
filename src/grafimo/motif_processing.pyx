# cython: profile=False, emit_code_comments=False, language_level=3

"""Functions to process the DNA motif Position Weight Matrix (PWM).

The PWM values are processed in order to obtain a scaled PWM and a 
corresponding P-value matrix, which is computed using a dynamic
programming algorithm (Staden, 1994)

The functions are written in Cython and for each one of them a Python
wrapper to call them is provided.

"""


from grafimo.utils import DNA_ALPHABET, lg2, RANGE, isListEqual, exception_handler
from grafimo.GRAFIMOException import NotValidAlphabetException, \
    FileReadError, NoDataFrameException, NotValidMotifMatrixException, \
    BGFileError, ScaledScoreMatrixException, MotifProcessingError
from grafimo.motif import Motif

from libc.stdlib cimport strtod

from typing import List, Dict, Set, Optional

import pandas as pd
import numpy as np

import os


#-------------------------- read the background file --------------------------#
cdef creadBGFile(bg_file, debug):
    """Read the background probabilties file. The background file must store the 
    probability for eachnucleotide (A, C, G, T).

    The background file must be given in Markov Background Model Format of 
    0th-order (http://meme-suite.org/doc/bfile-format.html).

    TODO: improve to read k-th order background files.

    ...

    Parameters
    ----------
    bg_file : str 
        path to the background file
    debug : bool
        trace the full error stack
        
    Returns
    -------
    Dict
        dictionary containing the background probabilities
    """

    if not isinstance(bg_file, str):
        errmsg = "Expected str, got %s.\n" % type(bg_file).__name__
        exception_handler(TypeError, errmsg, debug)
    if not os.path.isfile(bg_file):
        errmsg = "Unable to locate %s.\n" % bg_file
        exception_handler(FileNotFoundError, errmsg, debug)

    # ---- C vars
    cdef double prob  
    cdef char** endptr = NULL
    # ---- Py vars
    bg_dict: Dict = dict()
    found_nucs: Set = set()

    # bgfile parsing
    try:
        readlines = 0  # check if bgfile is empty
        ifstream = open(bg_file, mode="r") 
        for line in ifstream:
            if line[0] == "#": continue  # header
            if line[0].upper() in DNA_ALPHABET:
                nuc, prob_str = line.split()
                prob_str = prob_str.encode("UTF-8")  # from str to bytes
                prob = strtod(prob_str, endptr)
                assert prob > 0
                if nuc.upper() in found_nucs:
                    errmsg = "Found two times %s.\n" % nuc
                    exception_handler(BGFileError, errmsg, debug)
                bg_dict.update({nuc.upper(): prob})
                found_nucs.add(nuc.upper())
            else:
                errmsg = "Found symbol not part of the DNA alphabet: %s\n" % line[0]
                exception_handler(ValueError, errmsg, debug)
            # all nucs read
            if len(found_nucs) == len(DNA_ALPHABET): break
    except:
        errmsg = "An error occurred while parsing %s" % bg_file
        exception_handler(BGFileError, errmsg, debug)
    else:
        return bg_dict
    finally:
        ifstream.close()

# end creadBGfile()


def readBGfile(bg_file: str, debug: bool) -> Dict:
    """Python wrapper for creadBGfile() function 

    Read the background probabilties file. The background file must store the 
    probability for eachnucleotide (A, C, G, T).

    The background file must be given in Markov Background Model Format of 
    0th-order (http://meme-suite.org/doc/bfile-format.html).

    TODO: improve to read k-th order background files.

    ...

    Parameters
    ----------
    bg_file : str 
        path to the background file
    debug : bool
        trace the full error stack
        
    Returns
    -------
    Dict
        dictionary containing the background probabilities
    """

    return creadBGFile(bg_file, debug)


### uniform background distribution ###
cdef cget_uniformBG(alphabet, debug):
    """Compute a uniform background probability distribution for a given alphabet
    of characters.

    ...
        
    Parameters
    ----------
    alphabet : list 
        alphabet 
    debug : bool
        trace the full error stack
        
    Returns
    -------
    Dict
        dictionary containing a uniform probability background 
        distribution
    """

    if not isinstance(alphabet, list):
        errmsg = "Expected list, got %s.\n" % type(alphabet).__name__
        exception_handler(TypeError, errmsg, debug)
    if not isListEqual(alphabet, DNA_ALPHABET):
        errmsg = "The background alphabet is not the DNA alphabet.\n"
        exception_handler(ValueError, errmsg, debug)

    # ---- C vars
    cdef int alphasz  # alphabet size
    cdef double unifp  # uniform probabilities
    # ---- Py vars
    bg_dict: Dict = dict()

    alphasz = len(alphabet)
    unifp = 1.0 / <double>alphasz
    for i in range(alphasz): bg_dict.update({alphabet[i]: unifp})
    return bg_dict

# end cget_uniformBG()


def get_uniformBG(alphabet: List[str], debug: bool) -> Dict:
    """Python wrapper for cget_uniformBG() function 
    
    Compute a uniform background probability distribution for a given alphabet
    of characters.

    ...
        
    Parameters
    ----------
    alphabet : list 
        alphabet 
    debug : bool
        trace the full error stack
        
    Returns
    -------
    Dict
        dictionary containing a uniform probability background 
        distribution
    """

    return cget_uniformBG(alphabet, debug)


### post-process jaspar motif ###
cdef capply_pseudocount_jaspar(
    counts_matrix, 
    probs_matrix, 
    double pseudocount, 
    bgs,
    int width, 
    alphabet,
    nucsmap,
    debug
):
    """Apply pseudocount value to motif counts matrix. The motif PWM is given in
    JASPAR file format.

    ...

    Parameters
    ----------
    counts_matrix : numpy.ndarray
        motif raw counts matrix 
    probs_matrix : numpy.ndarray 
        motif probability motif matrix 
    pseudocount : numpy.double  
        pseudocount
    width : int
        motif width 
    bgs : dict 
        background probability distribution 
    alphabet : list 
        motif alphabet
    nucsmap : dict
        nucleotide ndarray index map
    debug : bool
        trace the full error stack
        
    Returns
    -------
    numpy.ndarray 
        processed motif probability matrix
    """

    if not isinstance(counts_matrix, np.ndarray):
        errmsg = "Expected numpy.ndarray, got %s.\n" % type(counts_matrix).__name__
        exception_handler(TypeError, errmsg, debug)
    if counts_matrix.size == 0 or sum(sum(counts_matrix)) == 0:
        errmsg = "Motif counts matrix is empty.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(probs_matrix, np.ndarray):
        errmsg = "Expected numpy.ndarray, got %s.\n" % type(probs_matrix).__name__
        exception_handler(TypeError, errmsg, debug)
    if probs_matrix.size == 0 or sum(sum(probs_matrix)) == 0:
        errmsg = "Motif probability matrix is empty.\n"
        exception_handler(ValueError, errmsg, probs_matrix)
    if not isinstance(alphabet, list):
        errmsg = "Expected list, got %s.\n" % type(alphabet).__name__
        exception_handler(TypeError, errmsg, debug)
    if not isListEqual(alphabet, DNA_ALPHABET):
        errmsg = "The motif is not built on DNA alphabet.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(bgs, dict):
        errmsg = "Expected dict, got %s.\n" % type(bgs).__name__
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(nucsmap, dict):
        errmsg = "Expected dict, got %s.\n" % type(nucsmap).__name__
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(pseudocount, float):
        errmsg = "Expected numpy.double, got %s.\n" % type(pseudocount).__name__
        exception_handler(TypeError, errmsg, debug)
    if pseudocount <= 0:
        errmsg = "Pseudocount values must be > 0.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(width, int):
        errmsg = "Expected int, got %s.\n" % type(width).__name__
        exception_handler(TypeError, errmsg, debug)
    if width <= 0:
        errmsg = "Forbidden motif width.\n"
        exception_handler(ValueError, errmsg, debug)

    # ---- C vars
    cdef double pseudo
    cdef int site_counts
    cdef double total_counts
    cdef double count
    cdef double bg

    pseudo = pseudocount
    proc_matrix = np.zeros(counts_matrix.shape, dtype=np.double)
    for j in range(width):
        site_counts = sum(counts_matrix[:, j])
        total_counts = <double>site_counts + pseudo
        for nuc in alphabet:
            nucidx = nucsmap[nuc]
            bg = bgs[nuc]
            assert bg > 0
            count = (
                (probs_matrix[nucidx, j] * <double>site_counts) + (pseudo * bg)
            )
            proc_matrix[nucidx, j] = count / total_counts
    assert proc_matrix.size != 0
    assert sum(sum(proc_matrix)) != 0
    return proc_matrix

# end capply_pseudocount_jaspar()


def apply_pseudocount_jaspar(
    counts_matrix: pd.DataFrame, 
    probs_matrix: pd.DataFrame, 
    pseudocount: np.double, 
    bgs: dict,
    width: int, 
    alphabet: List[str],
    nucsmap: Dict,
    debug: bool
) -> np.ndarray:
    """Python wrapper for capply_pseudocount_jaspar() function.
    
    Apply pseudocount value to motif counts matrix. The motif PWM is given in
    JASPAR file format.

    ...

    Parameters
    ----------
    counts_matrix : numpy.ndarray
        motif raw counts matrix 
    probs_matrix : numpy.ndarray 
        motif probability motif matrix 
    pseudocount : numpy.double  
        pseudocount
    width : int
        motif width 
    bgs : dict 
        background probability distribution 
    alphabet : list 
        motif alphabet
    nucsmap : dict
        nucleotide ndarray index map
    debug : bool
        trace the full error stack
        
    Returns
    -------
    numpy.ndarray 
        processed motif probability matrix
    """

    return capply_pseudocount_jaspar(
        counts_matrix, probs_matrix, pseudocount, bgs, width, alphabet, 
        nucsmap, debug
    )


### post-process meme motif ###
cdef capply_pseudocount_meme(
    probs_matrix, 
    double pseudocount, 
    int site_counts, 
    int width, 
    bgs, 
    alphabet,
    nucsmap,
    debug
):
    """Apply pseudocount value to motif counts matrix. The motif PWM is given in
    MEME file format.

    ...

    Parameters
    ----------
    probs_matrix : numpy.ndarray 
        motif probability motif matrix 
    pseudocount : float  
        pseudocount
    site_counts: int
        site counts
    width : int
        motif width 
    bgs : dict 
        background probability distribution 
    alphabet : list 
        motif alphabet
    nucsmap : dict
        nucleotide ndarray index map
    debug : bool
        trace the full error stack
        
    Returns
    -------
    numpy.ndarray 
        processed motif probability matrix
    """

    if not isinstance(probs_matrix, np.ndarray):
        errmsg = "Expected numpy.ndarray, got %s.\n" % type(probs_matrix).__name__
        exception_handler(TypeError, errmsg, debug)
    if probs_matrix.size == 0 or sum(sum(probs_matrix)) == 0:
        errmsg = "The probability matrix is empty.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(pseudocount, float):
        errmsg = "Expected float, got %s.\n" % type(pseudocount).__name__
        exception_handler(TypeError, errmsg, debug)
    if pseudocount <= 0:
        errmsg = "The pseudocount must be > 0."
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(site_counts, int):
        errmsg = "Expected int, got %s.\n" % type(site_counts).__name__
        exception_handler(TypeError, errmsg, debug)
    if site_counts <= 0:
        errmsg = "The site counts must be > 0.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(width, int):
        errmsg = "Expected int, got %s.\n" % type(width).__name__
        exception_handler(TypeError, errmsg, debug)
    if width <= 0:
        errmsg = "Forbidden motif width.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(bgs, dict):
        errmsg = "Excpected dict, got %s.\n" % type(bgs).__name__
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(alphabet, list):
        errmsg = "Expected list, got %s.\n" % type(alphabet).__name__
        exception_handler(TypeError, errmsg, debug)
    if not isListEqual(alphabet, DNA_ALPHABET):
        errmsg = "The motif is not built on DNA alphabet.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(nucsmap, dict):
        errmsg = "Expected dict, got %s.\n" % type(nucsmap).__name__
        exception_handler(TypeError, errmsg, debug)

    # ---- C vars
    cdef double total_counts
    cdef double bg
    cdef double count

    proc_matrix = np.zeros(probs_matrix.shape, dtype=np.double)
    total_counts = <double>site_counts + pseudocount
    for j in range(width):
        for nuc in alphabet:
            nucidx = nucsmap[nuc]
            bg = bgs[nuc]
            assert bg > 0
            count = ((probs_matrix[nucidx, j] * site_counts) + (pseudocount * bg))
            proc_matrix[nucidx, j] = count / total_counts

    return proc_matrix

# end of capply_pseudocount_meme()


def apply_pseudocount_meme(
    probs_matrix: np.ndarray, 
    pseudocount: np.double, 
    site_counts: int, 
    width: int, 
    bgs: Dict, 
    alphabet: List[str],
    nucsmap: Dict,
    debug: bool
) -> np.ndarray:
    """Python wrapper for capply_pseudocount_meme() function.
    
    Apply pseudocount value to motif counts matrix. The motif PWM is given in
    MEME file format.

    ...

    Parameters
    ----------
    probs_matrix : numpy.ndarray 
        motif probability motif matrix 
    pseudocount : float  
        pseudocount
    site_counts : int
        site counts
    width : int
        motif width 
    bgs : dict 
        background probability distribution 
    alphabet : list 
        motif alphabet
    nucsmap : dict
        nucleotide ndarray index map
    debug : bool
        trace the full error stack
        
    Returns
    -------
    numpy.ndarray 
        processed motif probability matrix
    """

    return capply_pseudocount_meme(
        probs_matrix, pseudocount, site_counts, width, bgs, alphabet, nucsmap, debug
    )


### compute log-odds ###
cdef ccompute_log_odds(probs_matrix, int width, bgs, alphabet, nucsmap, debug):
    """Computes motif log-odds scoring matrix from motif probability matrix.
    This same procedure is used for bith motifs PWMs in MEME or JASPAR formats.
        
    ...

    Parameters
    ----------
    probs_matrix : numpy.ndarray
        motif probability matrix 
    width : int 
        motif width
    bgs : dict 
        background probability distribution
    alphabet : list 
        motif alphabet 
    nucsmap : dict
        nucleotides index map
    debug : bool
        trace the full error stack
        
    Returns
    -------
    numpy.ndarray : 
        motif log-odds matrix
    """

    if not isinstance(probs_matrix, np.ndarray):
        errmsg = "Expected numpy.ndarray, got %s.\n" % type(probs_matrix).__name__
        exception_handler(TypeError, errmsg, debug)
    if probs_matrix.size == 0 or sum(sum(probs_matrix)) == 0:
        errmsg = "The motif probability matrix is empty.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(width, int):
        errmsg = "Expected int, got %s.\n" % type(width).__name__
        exception_handler(TypeError, errmsg, debug)
    if width <= 0:
        errmsg = "Forbidden motif width.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(bgs, dict):
        errmsg = "Expected dict, got %s.\n" % type(bgs).__name__
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(alphabet, list):
        errmsg = "Expected list, got %s.\n" % type(alphabet).__name__
        exception_handler(TypeError, errmsg, debug)
    if not isListEqual(alphabet, DNA_ALPHABET):
        errmsg = "The motif is not built on DNA alphabet.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(nucsmap, dict):
        errmsg = "Expected dict, got %s.\n" % type(nucsmap).__name__
        exception_handler(TypeError, errmsg, debug)

    # ---- C vars
    cdef double totBG
    cdef double totFG
    cdef double bg
    cdef double prob
    cdef double odds
    cdef double logodds
    cdef double epsilon = 0.001

    motif_log_odds = np.zeros(probs_matrix.shape, dtype=np.double)
    for nuc in alphabet:
        nucidx = nucsmap[nuc]
        bg = bgs[nuc]
        assert bg > 0
        totBG += bg
        for j in range(width):
            prob = probs_matrix[nucidx, j]
            assert prob > 0
            totFG += prob
            odds = prob / bg
            logodds = lg2(odds)
            motif_log_odds[nucidx, j] = logodds
    assert totBG - 1.0 < epsilon
    assert totFG - width < epsilon

    return motif_log_odds

# end of ccompute_log_odds()


def compute_log_odds(
    probs_matrix: pd.DataFrame, 
    width: int, 
    bgs: Dict, 
    alphabet: List[str],
    nucsmap: dict,
    debug: bool
) -> np.ndarray:
    """Python wrapper function for compute_lo_odds() function 
    
    Computes motif log-odds scoring matrix from motif probability matrix.
    This same procedure is used for bith motifs PWMs in MEME or JASPAR formats.
        
    ...

    Parameters
    ----------
    probs_matrix : numpy.ndarray
        motif probability matrix 
    width : int 
        motif width
    bgs : dict 
        background probability distribution
    alphabet : list 
        motif alphabet 
    nucsmap : dict
        nucleotides index map
    debug : bool
        trace the full error stack
        
    Returns
    -------
    numpy.ndarray : 
        motif log-odds matrix
    """

    return ccompute_log_odds(probs_matrix, width, bgs, alphabet, nucsmap, debug)


### DP p-value matrix computation ###
cdef ccomp_pval_mat(motif, debug):
    """Computes the motif P-value matrix using the dynamic programming 
    algorithm presented in Staden R, 1994.
    
    The P-value matrix allows to compute in constant time the statistical 
    significance of the motif occurrence candidates log-likelihood scores.

    ...
        
    Parameters
    ----------
    motif : Motif 
        motif object
    debug : bool
        trace the full error stack
        
    Returns
    -------
    numpy.ndarray
        P-value matrix
    """

    if not isinstance(motif, Motif):
        errmsg = "Expected Motif, got %s.\n" % type(motif).__name__
        exception_handler(TypeError, errmsg, debug)
    if not motif.isScaled:
        errmsg = "The motif score matrix has not been scaled yet.\n"
        exception_handler(MotifProcessingError, errmsg, debug)

    # ---- C vars
    cdef int width  # motif width
    cdef double bg  # background 
    cdef double source  

    # ---- Py vars
    score_matrix = motif.scoreMatrix
    width = motif.width
    alphabet = motif.alphabet
    bgs = motif.bg
    nucsmap = motif.nucsmap

    pval_mat = np.zeros((width, (RANGE * width + 1)))
    for pos in range(width):
        for nuc in alphabet:
            bg = bgs[nuc]
            assert bg > 0
            if pos == 0:
                pval_mat[0, score_matrix[nucsmap[nuc], pos]] += np.double(1 * bg)
            else:
                idxs = np.where(pval_mat[pos-1, :] > 0)[0]
                for idx in idxs:
                    if pval_mat[pos-1, idx] != 0:
                        source = pval_mat[pos-1, idx]
                        pval_mat[pos, score_matrix[nucsmap[nuc], pos] + idx] += source * bg
    # keep only last row
    pval_mat = pval_mat[width-1]

    return pval_mat

# end of ccomp_pval_mat()


def comp_pval_mat(motif: Motif, debug: bool) -> np.ndarray:
    """Python wrapper for comp_pval_mat() function 
    
    Computes the motif P-value matrix using the dynamic programming 
    algorithm presented in Staden R, 1994.
    
    The P-value matrix allows to compute in constant time the statistical 
    significance of the motif occurrence candidates log-likelihood scores.

    ...
        
    Parameters
    ----------
    motif : Motif 
        motif object
    debug : bool
        trace the full error stack
        
    Returns
    -------
    numpy.ndarray
        P-value matrix
    """

    return ccomp_pval_mat(motif, debug)

