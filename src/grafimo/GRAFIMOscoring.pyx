# cython: profile=False, emit_code_comments=False, language_level=3

"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Script that contains Cython utilities functions for the
k-mer scoring step

"""

import pandas as pd
from grafimo.GRAFIMOException import ValueException, NoDataFrameException
from grafimo.utils import die, REV_COMPL

### score each k-mer ###
cdef cscore_seq(seq, scoreMatrix, pval_mat, minScore, scale, width, offset):
    """
        Score the given sequence with values in the given score matrix
        ----
        Params: 
            seq (str) : sequence to score
            scoreMatrix (pd.DataFrame) : score matrix
            pval_mat (np.array) : array of the p-values
            minScore (int) : minimum score of the score matrix
            scale (int) : scale factor used 
            width (int) : motif's width 
            offset (double) : offset of the scaled score
        ----
        Returns:
            logodds_score (double) : log-odds score of the sequence
            pvalue (double) : p-value of the obtained score
    """

    if not isinstance(seq, str):
        raise ValueException("Cannot score the sequence")
        die(1)

    if not isinstance(scoreMatrix, pd.DataFrame):
        raise NoDataFrameException("Cannot score sequence with the given score matrix")

    assert minScore != None
    assert width > 0

    cdef int score = 0
    cdef int rev_score = 0
    cdef double pvalue
    cdef double logodds_score

    for idx, nuc in enumerate(seq):
        nuc = nuc.upper()

        if nuc == 'N':
            score = minScore
            break # we don't go further

        score += scoreMatrix.loc[nuc, idx]

    assert score >= 0 and score <= width*1000 # score must be in range
    
    pvalue = getPvalue(pval_mat, score)

    # turns back score to log-odds
    logodds_score = recover_logodds(score, scale, width, offset)

    return logodds_score, pvalue

def score_seq(seq, scoreMatrix, pval_mat, minScore, scale, width, offset):
    """ Python call function """

    return cscore_seq(seq, scoreMatrix, pval_mat, minScore, scale, width, offset)

### get the p-value ###
cdef cgetPvalue(pval_mat, score):
    """
        Compute the p-value of a given score
        ----
        Params:
            pval_mat (np.array) : array of p-values
            score (int) : input score
        ----
        Returns:
            pvalue (double) : p-value for the input score
    """

    if not isinstance(score, int):
        raise ValueException("The score must be an int. Cannot compute p-value")
        die(1)

    assert sum(pval_mat) > 0
    assert score >= 0

    cdef double tot
    cdef double pvalue

    tot=sum(pval_mat[:])
    pvalue=(sum(pval_mat[score:]))/tot

    return pvalue

def getPvalue(pval_mat, score):
    """ Python call function """

    return cgetPvalue(pval_mat, score)

### recover log-odds score ###
cdef crecover_logodds(score, scale, width, offset):
    """
        Recover the log-odds score, given the scaled one
        ----
        Params: 
            score (int) : scaled score 
            scale (int) : scaling factor used
            width (int) : motif width
            offset (double) : offset of the scaled score
        ----
        Returns:
             logodds (double) : log-odds representation of the scaled score
    """

    assert score >= 0
    assert scale > 0

    cdef double logodds

    logodds=(score/scale)+(width*offset)

    return logodds

def recover_logodds(score, scale, width, offset):
    """ Python call function """

    return crecover_logodds(score, scale, width, offset)
