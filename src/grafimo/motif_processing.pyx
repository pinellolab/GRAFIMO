# cython: profile=False, emit_code_comments=False, language_level=3

"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Cython script that contains some functions to process the motif
matrix

"""


from libc.stdlib cimport strtod
from grafimo.utils import DNA_ALPHABET, lg2, RANGE, isListEqual
from grafimo.GRAFIMOException import NotValidAlphabetException, FileReadingException, NoDataFrameException, \
                                        NotValidMotifMatrixException, NotValidBGException, ScaledScoreMatrixException
import pandas as pd
import numpy as np


### read the background file ###
cdef creadBGFile(bg_file):
    """
        Read the background file given by the user
        ----
        Parameters:
            bg_file (str) : path to the background file
        ----
        Returns:
            bg_dict (dict) : dictionary that contains the background probabilities
                             defined in the input file
    """

    cdef double prob
    cdef char** endptr = NULL

    bg_dict = {}
    found_nucs = set()

    try:
        bgf = open(bg_file, mode='r') # open the file in read only mode

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

            if len(found_nucs) == len(DNA_ALPHABET): # read all the nucleotides
                if found_nucs != set(DNA_ALPHABET): # some symbol was read twice or more
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


def readBGfile(bg_file):
    """ Python call function """

    return creadBGFile(bg_file)


### uniform background distribution ###
cdef cget_uniformBG(alphabet):
    """
        Returns a uniform probability distribution for a given alphabet
        ----
        Parameters:
            alphabet (list) : alphabet in list form
        ----
            Returns:
                bg_dict (dict) : dictionary of the uniform probability
    """

    if not isinstance(alphabet, list):
        raise NotValidAlphabetException("\n\nERROR: the given alphabet is not in a list")

    if not isListEqual(alphabet, DNA_ALPHABET):
        raise NotValidAlphabetException("\n\nERROR: the given alphabet is not a valid DNA alphabet")

    cdef int alpha_len # length of the alphabet
    cdef double unifp # uniform probability value

    bg_dict = {}

    alpha_len = len(alphabet)
    unifp = 1.0 / <double>alpha_len

    for i in range(alpha_len):
        bg_dict.update({alphabet[i]: unifp})

    return bg_dict
# end cget_uniformBG()


def get_uniformBG(alphabet):
    """ Python call function """

    return cget_uniformBG(alphabet)


### post-process jaspar motif ###
cdef capply_pseudocount_jaspar(counts_matrix, probs_matrix, double pseudocount, bgs,
                                    int width, alphabet):
    """
        Apply the pseudocount to raw counts matrix (JASPAR motif file)
        ----
        Params:
            counts_matrix (pd.DataFrame) : raw counts matrix 
            probs_matrix (pd.DataFrame) : probability matrix 
            pseudocount (np.double) : pseudocount to apply to the motif
            width (int) : motif width 
            bgs (dict) : background distribution
            alphabet (list) : alphabet of the motif
        ----
        Returns:
            proc_matrix (pd.DataFrame) : resulting probability matrix 
                                            (pseudocount applied)
    """

    if not isinstance(counts_matrix, pd.DataFrame):
        errmsg = "\n\nERROR: the given motif matrix must be an instance of pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if counts_matrix.empty:
        raise NotValidMotifMatrixException("\n\nERROR: the given motif matrix is empty")

    if not isinstance(probs_matrix, pd.DataFrame):
        errmsg = "\n\nERROR: the given motif matrix must be an instance of pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if probs_matrix.empty:
        raise NotValidMotifMatrixException("\n\nERROR: the given motif matrix is empty")

    if not isinstance(alphabet, list):
        raise NotValidAlphabetException("\n\nERROR: the given alphabet is not in a list")

    if not isListEqual(alphabet, DNA_ALPHABET):
        raise NotValidAlphabetException("\n\nERROR: the given alphabet is not a valid DNA alphabet")

    if not isinstance(bgs, dict):
        raise NotValidBGException("\n\nERROR: the background must be in a dict data-structure")

    assert pseudocount > 0
    assert width > 0

    cdef double pseudo
    cdef int site_counts
    cdef double total_counts
    cdef double count
    cdef double bg

    pseudo = pseudocount

    proc_matrix = pd.DataFrame(index=list(counts_matrix.index), columns=list(counts_matrix.columns),
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


def apply_pseudocount_jaspar(counts_matrix, probs_matrix, pseudocount, bgs,
                                width, alphabet):
    """ Python call function """

    return capply_pseudocount_jaspar(counts_matrix, probs_matrix, pseudocount, bgs,
                                        width, alphabet)


### post-process meme motif ###
cdef capply_pseudocount_meme(probs_matrix, double pseudocount, int site_counts, int width, bgs, alphabet):
    """
        Add pseudocount to the motif matrix (MEME format)
        ----
        Params:
            probs_matrix (pd.DataFrame) : probability matrix of the motif
            pseudocount (np.double) : pseudocount to apply to the motif matrix
            site_counts (int) : site counts of the motif
            width (int) : motif width
            bgs (dict) : background distribution 
            alphabet (list) : motif alphabet
        ----
        Returns:
             proc_motif (pd.DataFrame) : processed probability matrix
    """

    if not isinstance(probs_matrix, pd.DataFrame):
        errmsg = "\n\nERROR: the given motif matrix must be an instance of pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if probs_matrix.empty:
        raise NotValidMotifMatrixException("\n\nERROR: the given motif matrix is empty")

    if not isinstance(alphabet, list):
        raise NotValidAlphabetException("\n\nERROR: the given alphabet is not in a list")

    if not isListEqual(alphabet, DNA_ALPHABET):
        raise NotValidAlphabetException("\n\nERROR: the given alphabet is not a valid DNA alphabet")

    if not isinstance(bgs, dict):
        raise NotValidBGException("\n\nERROR: the background must be in a dict data-structure")

    assert pseudocount > 0
    assert site_counts > 0
    assert width > 0

    cdef double total_counts
    cdef double bg
    cdef double count

    proc_matrix = pd.DataFrame(index=list(probs_matrix.index), columns=list(probs_matrix.columns),
                                data=np.double(0))

    total_counts = <double>site_counts + pseudocount

    for j in range(width):
        for nuc in alphabet:

            bg = bgs[nuc]
            count = ((probs_matrix.loc[nuc, j]*site_counts)+
                     (pseudocount*bg))
            proc_matrix.loc[nuc, j] = count / total_counts
        # end for
    # end for

    return proc_matrix
# end of capply_pseudocount_meme()


def apply_pseudocount_meme(probs_matrix, pseudocount, site_counts, width, bgs, alphabet):
    """ Python call function """

    return capply_pseudocount_meme(probs_matrix, pseudocount, site_counts, width,
                                   bgs, alphabet)


### compute log-odds ###
cdef ccompute_log_odds(probs_matrix, int width, bgs, alphabet):
    """
        Computes the log-odds scores for the motif probability matrix
        ----
        Params:
            probs_matrix (pd.DataFrame) : motif probability matrix 
            width (int) : width of the motif
            bgs (dict) : background distribution
            alphabet (list) : alphabet of the motif 
        ----
        Returns:
             motif_log_odds (pd.DataFrame) : log-odds motif matrix
    """

    if not isinstance(probs_matrix, pd.DataFrame):
        errmsg = "\n\nERROR: the given motif matrix must be an instance of pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if probs_matrix.empty:
        raise NotValidMotifMatrixException("\n\nERROR: the given motif matrix is empty")

    if not isinstance(alphabet, list):
        raise NotValidAlphabetException("\n\nERROR: the given alphabet is not in a list")

    if not isinstance(bgs, dict):
        raise NotValidBGException("\n\nERROR: the background must be in a dict data-structure")

    if not isListEqual(alphabet, DNA_ALPHABET):
        raise NotValidAlphabetException("\n\nERROR: the given alphabet is not a valid DNA alphabet")

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


def compute_log_odds(probs_matrix, width, bgs, alphabet):
    """ Python call function """

    return ccompute_log_odds(probs_matrix, width, bgs, alphabet)


### DP p-value matrix computation ###
cdef ccomp_pval_mat(motif):
    """
        Computes the p-value matrix needed for the computation of the p-value,
        using a Dynamic-Programming algorithm
        ----
        Parameters:
            motif (Motif) : motif object
        ----
        Returns:
            pval_mat (np.array) : p-value matrix computed from the motif scores
    """

    if not motif.getIsScaled():
        errmsg = "\n\nERROR: the motif has no scaled score matrix. Cannot compute the p-value matrix"
        raise ScaledScoreMatrixException(errmsg)

    if not isinstance(motif.getMotif_scoreMatrix(), pd.DataFrame):
        errmsg = "\n\nERROR: the given motif matrix must be an instance of pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if motif.getMotif_scoreMatrix().empty:
        raise NotValidMotifMatrixException("\n\nERROR: the given motif matrix is empty")

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
                idxs = np.where(pval_mat[pos - 1, :] > 0)[0]

                for idx in idxs:
                    if pval_mat[pos-1, idx] != 0:
                        source = pval_mat[pos-1, idx]
                        pval_mat[pos, score_matrix.loc[nuc, pos] + idx] += source*bg
                    # end if
                # end for
            # end if
        # end for
    # end for

    # get the interesting part of the matrix
    pval_mat = pval_mat[width-1]

    return pval_mat
# end of ccomp_pval_mat()


def comp_pval_mat(motif):
    """ Python call function """

    return ccomp_pval_mat(motif)

