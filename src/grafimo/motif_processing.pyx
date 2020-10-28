# cython: profile=False, emit_code_comments=False, language_level=3

"""Functions to process the DNA motif Position Weight Matrix (PWM).

The PWM values are processed in order to obtain a scaled PWM and a 
corresponding P-value matrix, which is computed using a dynamic
programming algorithm (Staden, 1994)

The functions are written in Cython and for each one of them a Python
wrapper to call them is provided.

"""


from grafimo.utils import DNA_ALPHABET, lg2, RANGE, isListEqual
from grafimo.GRAFIMOException import NotValidAlphabetException, \
    FileReadingException, NoDataFrameException, NotValidMotifMatrixException, \
    NotValidBGException, ScaledScoreMatrixException
from grafimo.motif import Motif
from typing import List, Dict, Optional
from libc.stdlib cimport strtod
import pandas as pd
import numpy as np


### read the background file ###
cdef creadBGFile(bg_file):
    """Read the background file, which contains the background 
    frequencies for each nucleotide.

    The background file is given in Markov Background Model Format
    (http://meme-suite.org/doc/bfile-format.html).

    Parameters
    ----------
    bg_file : str 
        path to the background file
        
    Returns
    -------
    Dict
        dictionary containing the background probabilities
    """

    cdef double prob
    cdef char** endptr = NULL

    bg_dict: dict
    found_nucs: set
    errmsg: str

    bg_dict = dict()
    found_nucs = set()

    try:
        bgf = open(bg_file, mode='r') 

        for line in bgf:

            if line[0] in DNA_ALPHABET:
                nuc, prob_str = line.split('\t')
                prob_str = prob_str.encode("UTF-8") # from str to bytes
                prob = strtod(prob_str, endptr)

                assert prob > 0

                bg_dict.update({nuc: prob})
                found_nucs.add(nuc)

            else:
                errmsg = "\n\nERROR: the read symbol is not part of the DNA alphabet"
                raise NotValidAlphabetException(errmsg)

            if len(found_nucs) == len(DNA_ALPHABET):
                # read all the nucleotides 
                if found_nucs != set(DNA_ALPHABET): 
                    # some symbol was read twice or more
                    errmsg = "\n\nERROR: some alphabet symbols were read twice or more"
                    raise NotValidAlphabetException(errmsg)

                break
            # end if
        # end for

    except:
        errmsg=' '.join(["Unable to read file", bg_file])
        raise FileReadingException(errmsg)

    else:
        return bg_dict

    finally:
        bgf.close() # close the file
    # end try
# end creadBGfile()


def readBGfile(bg_file: str) -> Dict:
    """Python wrapper for creadBGfile() function 

    Read the background file, which contains the background 
    frequencies for each nucleotide.

    The background file is given in Markov Background Model Format
    (http://meme-suite.org/doc/bfile-format.html).

    Parameters
    ----------
    bg_file : str 
        path to the background file
        
    Returns
    -------
    Dict
        dictionary containing the background probabilities
    """

    return creadBGFile(bg_file)


### uniform background distribution ###
cdef cget_uniformBG(alphabet):
    """Compute a uniform probability background distribution for a given 
    alphabet of characters
        
    Parameters
    ----------
    alphabet : list 
        alphabet 
        
    Returns
    -------
    Dict
        dictionary containing a uniform probability background 
        distribution
    """

    errmsg: str
    if not isinstance(alphabet, list):
        errmsg = "\n\nERROR: the given alphabet is not in a list"
        raise NotValidAlphabetException(errmsg)

    if not isListEqual(alphabet, DNA_ALPHABET):
        errmsg = "\n\nERROR: the given alphabet is not a valid DNA alphabet"
        raise NotValidAlphabetException(errmsg)

    cdef int alpha_len # length of the alphabet
    cdef double unifp # uniform probability value

    bg_dict: Dict

    bg_dict = dict()

    alpha_len = len(alphabet)
    unifp = 1.0 / <double>alpha_len

    for i in range(alpha_len):
        bg_dict.update({alphabet[i]: unifp})

    return bg_dict
# end cget_uniformBG()


def get_uniformBG(alphabet: List[str]) -> Dict:
    """Python wrapper for cget_uniformBG() function 
    
    Compute a uniform probability background distribution for a given 
    alphabet of characters
        
    Parameters
    ----------
    alphabet : list 
        alphabet 
        
    Returns
    -------
    Dict
        dictionary containing a uniform probability background 
        distribution
    """

    return cget_uniformBG(alphabet)


### post-process jaspar motif ###
cdef capply_pseudocount_jaspar(
    counts_matrix, 
    probs_matrix, 
    double pseudocount, 
    bgs,
    int width, 
    alphabet
):
    """Apply the pseudocount value to the raw counts matrix PWM for the 
    motif given in JASPAR file format.

    Parameters
    ----------
    counts_matrix : pandas.DataFrame
        raw counts PWM matrix 
    probs_matrix : pandas.DataFrame 
        probability PWM matrix 
    pseudocount : numpy.double  
        pseudocount to add to the raw counts PWM motif
    width : int
        motif width 
    bgs : dict 
        background probability distribution of the DNA alphabet
    alphabet : list 
        DNA motif alphabet
        
    Returns
    -------
    pandas.DataFrame 
        resulting probability matrix with the pseudocount applied
    """

    errmsg: str
    if not isinstance(counts_matrix, pd.DataFrame):
        errmsg = "\n\nERROR: the given motif matrix must be an instance of pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if counts_matrix.empty:
        errmsg = "\n\nERROR: the given motif matrix is empty"
        raise NotValidMotifMatrixException(errmsg)

    if not isinstance(probs_matrix, pd.DataFrame):
        errmsg = "\n\nERROR: the given motif matrix must be an instance of pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if probs_matrix.empty:
        errmsg = "\n\nERROR: the given motif matrix is empty"
        raise NotValidMotifMatrixException(errmsg)

    if not isinstance(alphabet, list):
        errmsg = "\n\nERROR: the given alphabet is not in a list"
        raise NotValidAlphabetException(errmsg)

    if not isListEqual(alphabet, DNA_ALPHABET):
        errmsg = "\n\nERROR: the given alphabet is not a valid DNA alphabet"
        raise NotValidAlphabetException(errmsg)

    if not isinstance(bgs, dict):
        errmsg = "\n\nERROR: the background must be in a dict data-structure"
        raise NotValidBGException(errmsg)

    assert pseudocount > 0
    assert width > 0

    cdef double pseudo
    cdef int site_counts
    cdef double total_counts
    cdef double count
    cdef double bg

    pseudo = pseudocount

    proc_matrix = pd.DataFrame(index=list(counts_matrix.index), 
                               columns=list(counts_matrix.columns), 
                               data=np.double(0))

    for j in range(width):

        site_counts = sum(counts_matrix.loc[:, j])
        total_counts = <double>site_counts+pseudo

        for nuc in alphabet:

            bg = bgs[nuc]

            assert bg > 0

            count = ((probs_matrix.loc[nuc, j] * <double>site_counts) +
                        (pseudo*bg))
            proc_matrix.loc[nuc, j] = count / total_counts
        # end for
    # end for

    return proc_matrix
# end capply_pseudocount_jaspar()


def apply_pseudocount_jaspar(
    counts_matrix: pd.DataFrame, 
    probs_matrix: pd.DataFrame, 
    pseudocount: np.double, 
    bgs: dict,
    width: int, 
    alphabet: List[str]
) -> pd.DataFrame:
    """Python wrapper for capply_pseudocount_jaspar() function 
    
    Apply the pseudocount value to the raw counts matrix PWM for the 
    motif given in JASPAR file format.

    Parameters
    ----------
    counts_matrix : pandas.DataFrame
        raw counts PWM matrix 
    probs_matrix : pandas.DataFrame 
        probability PWM matrix 
    pseudocount : numpy.double  
        pseudocount to add to the raw counts PWM motif
    width : int
        motif width 
    bgs : dict 
        background probability distribution of the DNA alphabet
    alphabet : list 
        DNA motif alphabet
        
    Returns
    -------
    pandas.DataFrame 
        resulting probability matrix with the pseudocount applied
    """

    return capply_pseudocount_jaspar(counts_matrix, probs_matrix, pseudocount, 
                                     bgs, width, alphabet)


### post-process meme motif ###
cdef capply_pseudocount_meme(
    probs_matrix, 
    double pseudocount, 
    int site_counts, 
    int width, 
    bgs, 
    alphabet
):
    """Add pseudocount value to the probability motif PWM matrix in MEME 
    file format.
        
    Parameters:
    probs_matrix : pandas.DataFrame
        motif probability PWM matrix
    pseudocount : numpy.double
        pseudocount to apply to the motif PWM matrix
    site_counts : int 
        motif site counts motif
    width : int 
        motif width
    bgs : dict 
        background probability distribution for the DNA motif alphabet 
    alphabet : list 
        DNA motif alphabet
        
    Returns
    -------
    pandas.DataFrame 
        resulting probability matrix with the pseudocount applied
    """

    errmsg: str
    if not isinstance(probs_matrix, pd.DataFrame):
        errmsg = "\n\nERROR: the given motif matrix must be an instance of "
        errmsg += "pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if probs_matrix.empty:
        errmsg = "\n\nERROR: the given motif matrix is empty"
        raise NotValidMotifMatrixException(errmsg)

    if not isinstance(alphabet, list):
        errmsg = "\n\nERROR: the given alphabet is not in a list"
        raise NotValidAlphabetException(errmsg)

    if not isListEqual(alphabet, DNA_ALPHABET):
        errmsg = "\n\nERROR: the given alphabet is not a valid DNA alphabet"
        raise NotValidAlphabetException(errmsg)

    if not isinstance(bgs, dict):
        errmsg = "\n\nERROR: the background must be in a dict data-structure"
        raise NotValidBGException(errmsg)

    assert pseudocount > 0
    assert site_counts > 0
    assert width > 0

    cdef double total_counts
    cdef double bg
    cdef double count

    proc_matrix = pd.DataFrame(index=list(probs_matrix.index), 
                               columns=list(probs_matrix.columns),  
                               data=np.double(0))

    total_counts = <double>site_counts + pseudocount

    for j in range(width):
        for nuc in alphabet:

            bg = bgs[nuc]
            count = ((probs_matrix.loc[nuc, j] * site_counts) + 
                        (pseudocount * bg))
            proc_matrix.loc[nuc, j] = count / total_counts
        # end for
    # end for

    return proc_matrix
# end of capply_pseudocount_meme()


def apply_pseudocount_meme(
    probs_matrix: pd.DataFrame, 
    pseudocount: np.double, 
    site_counts: int, 
    width: int, 
    bgs: Dict, 
    alphabet: List[str]
) -> pd.DataFrame:
    """Python wrapper for capply_pseudocount_meme() function 
    
    Add pseudocount value to the probability motif PWM matrix in MEME 
    file format.
        
    Parameters:
    probs_matrix : pandas.DataFrame
        motif probability PWM matrix
    pseudocount : numpy.double
        pseudocount to apply to the motif PWM matrix
    site_counts : int 
        motif site counts motif
    width : int 
        motif width
    bgs : dict 
        background probability distribution for the DNA motif alphabet 
    alphabet : list 
        DNA motif alphabet
        
    Returns
    -------
    pandas.DataFrame 
        resulting probability matrix with the pseudocount applied
    """

    return capply_pseudocount_meme(probs_matrix, pseudocount, site_counts, 
                                   width, bgs, alphabet)


### compute log-odds ###
cdef ccompute_log_odds(probs_matrix, int width, bgs, alphabet):
    """Computes the log-odds scores from the values of the motif  
    probability matrix.

    This step is madatory either if the motif PWM has been given in
    JASPAR or MEME file format 
        
    Parameters
    ----------
    probs_matrix : pd.DataFrame
        motif probability matrix 
    width : int 
        motif width
    bgs : dict 
        background probability distribution
    alphabet : list 
        DNA motif alphabet 
        
    Returns
    -------
    pandas.DataFrame : 
        motif matrix with log-odds values
    """

    errmsg: str
    if not isinstance(probs_matrix, pd.DataFrame):
        errmsg = "\n\nERROR: the given motif matrix must be an instance of pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if probs_matrix.empty:
        errmsg = "\n\nERROR: the given motif matrix is empty"
        raise NotValidMotifMatrixException(errmsg)

    if not isinstance(alphabet, list):
        errmsg = "\n\nERROR: the given alphabet is not in a list"
        raise NotValidAlphabetException(errmsg)

    if not isinstance(bgs, dict):
        errmsg = "\n\nERROR: the background must be in a dict data-structure"
        raise NotValidBGException(errmsg)

    if not isListEqual(alphabet, DNA_ALPHABET):
        errmsg = "\n\nERROR: the given alphabet is not a valid DNA alphabet"
        raise NotValidAlphabetException(errmsg)

    assert width > 0

    cdef double totBG
    cdef double totFG
    cdef double bg
    cdef double prob
    cdef double odds
    cdef double logodds
    cdef double epsilon = 0.001

    motif_log_odds = pd.DataFrame(index=list(probs_matrix.index), columns=list(probs_matrix.columns),
                                    data=np.double(0))

    for nuc in alphabet:

        bg = bgs[nuc]

        assert bg > 0

        totBG += bg

        for j in range(width):

            prob = probs_matrix.loc[nuc, j]

            assert prob > 0

            totFG += prob
            odds = prob / bg
            logodds = lg2(odds)

            motif_log_odds.loc[nuc, j] = logodds
        # end for
    # end for

    assert totBG - 1.0 < epsilon
    assert totFG - width < epsilon

    return motif_log_odds
# end of ccompute_log_odds()


def compute_log_odds(
    probs_matrix: pd.DataFrame, 
    width: int, 
    bgs: Dict, 
    alphabet: List[str]
) -> pd.DataFrame:
    """Python wrapper function for compute_lo_odds() function 
    
    Computes the log-odds scores from the values of the motif  
    probability matrix.

    This step is madatory either if the motif PWM has been given in
    JASPAR or MEME file format 
        
    Parameters
    ----------
    probs_matrix : pd.DataFrame
        motif probability matrix 
    width : int 
        motif width
    bgs : dict 
        background probability distribution
    alphabet : list 
        DNA motif alphabet 
        
    Returns
    -------
    pandas.DataFrame : 
        motif matrix with log-odds values
    """

    return ccompute_log_odds(probs_matrix, width, bgs, alphabet)


### DP p-value matrix computation ###
cdef ccomp_pval_mat(motif):
    """Computes the P-value matrix using the dynamic programming 
    algorithm presented in Staden, 1994.
    
    The P-value matrix is required to assign a P-value in constant time 
    to each log-likelihood score for the motif occurrence candidates
        
    Parameters
    ----------
    motif : Motif 
        motif object
        
    Returns
    -------
    numpy.array
        P-value matrix
    """

    errmsg: str
    if not motif.getIsScaled():
        errmsg = "\n\nERROR: the motif has no scaled score matrix. Cannot compute the p-value matrix"
        raise ScaledScoreMatrixException(errmsg)

    if not isinstance(motif.getMotif_scoreMatrix(), pd.DataFrame):
        errmsg = "\n\nERROR: the given motif matrix must be an instance of pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if motif.getMotif_scoreMatrix().empty:
        errmsg = "\n\nERROR: the given motif matrix is empty"
        raise NotValidMotifMatrixException(errmsg)

    cdef int width
    cdef double bg
    cdef double source

    score_matrix = motif.getMotif_scoreMatrix()

    width = motif.getWidth()
    alphabet = motif.getAlphabet()
    bgs = motif.getBg()

    pval_mat = np.zeros((width, RANGE * width + 1))

    for pos in range(width):
        for nuc in alphabet:

            bg = bgs[nuc]

            if pos == 0:
                pval_mat[0, score_matrix.loc[nuc, pos]] += np.double(1 * bg)
            else:
                idxs = np.where(pval_mat[pos-1, :] > 0)[0]

                for idx in idxs:
                    if pval_mat[pos-1, idx] != 0:
                        source = pval_mat[pos-1, idx]
                        pval_mat[pos, score_matrix.loc[nuc, pos] + idx] += source * bg
                    # end if
                # end for
            # end if
        # end for
    # end for

    # get the interesting part of the matrix
    pval_mat = pval_mat[width-1]

    return pval_mat
# end of ccomp_pval_mat()


def comp_pval_mat(motif: Motif) -> np.array:
    """Python wrapper for comp_pval_mat() function 
    
    Computes the P-value matrix using the dynamic programming 
    algorithm presented in Staden, 1994.
    
    The P-value matrix is required to assign a P-value in constant time 
    to each log-likelihood score for the motif occurrence candidates
        
    Parameters
    ----------
    motif : Motif 
        motif object
        
    Returns
    -------
    numpy.array
        P-value matrix
    """

    return ccomp_pval_mat(motif)

