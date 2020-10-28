"""Definition of ResultTmp class.

The class is used to store intermediate results during the motif 
occurrences candidates scoring step of GRAFIMO analysis.
"""

from grafimo.GRAFIMOException import ValueException 
from typing import List
import numpy as np


class ResultTmp(object):
    """
    This class stores intermediate results of the scoring step of 
    GRAFIMO analysis.

    The fields of a ResultTmp instance resemble those of the final
    TSV and HTML reports.

    A ResultTmp instance is created for each thread created when scoring 
    the motif occurrences. The different instances are then merged in
    a single results report, which will be printed in the three output
    files or directly on terminal.

    ...

    Attributes
    ----------
    _seqnames : list
        sequence names
    _seqs : list
        sequences
    _chroms : list
        chromosomes
    _starts : list
        sequence start coordinates
    _stops : list
        sequence stop coordinates
    _strands : list
        sequence DNA strand
    _scores : list
        sequence scores
    _pvalues : list
        sequence P-values
    _frequencies : list
        sequence frequencies
    _references : list
        flag value to state if the sequences contain genomic variants
    
    Methods
    -------
    get_seqnames()
        return the sequence names
    get_seqs()
        return the sequences
    get_chroms()
        return the sqeuence chromosomes
    get_starts()
        return the sequence start coordinates
    get_stops()
        return the sequence stop coordinates
    get_strands()
        return the sequence DNA strands
    get_scores()
        return the sequence scores
    get_pvalues()
        return the sequence P-values
    get_frequencies()
        return the sequence frequencies
    get_references()
        return flag values for sequences, stating if they contain or not
        genomic variants
    """

    #-------------------------------------------------------------------
    # ResultTmp attributes
    #-------------------------------------------------------------------
    _seqnames: List[str]
    _seqs: List[str]
    _chroms: List[str]
    _starts: List[int]
    _stops: List[int]
    _strands: List[str]
    _scores: List[np.double]
    _pvalues: List[np.double]
    _frequencies: List[int]
    _references: List[str]


    #-------------------------------------------------------------------
    # ResultTmp methods
    #-------------------------------------------------------------------
    def __init__(self,
                 seqnames: List[str],
                 seqs: List[str],
                 chroms: List[str],
                 starts: List[int],
                 stops: List[int],
                 strands: List[str],
                 scores: List[np.double],
                 pvalues: List[np.double],
                 frequencies: List[int],
                 references: List[str]
    ):

        assert len(seqnames) == len(seqs)
        assert len(seqnames) == len(chroms)
        assert len(seqnames) == len(starts)
        assert len(seqnames) == len(stops)
        assert len(seqnames) == len(strands)
        assert len(seqnames) == len(scores)
        assert len(seqnames) == len(pvalues)
        assert len(seqnames) == len(frequencies)
        assert len(seqnames) == len(references)

        errmsg: str 
        
        errmsg = ' '.join(["\n\nERROR: unable to store temporary results.",
                                "Wrong data-type given"]) 
        if not isinstance(seqnames, list):
            raise ValueError(errmsg)

        if not isinstance(seqs, list):
            raise ValueError(errmsg)

        if not isinstance(chroms, list):
            raise ValueError(errmsg)

        if not isinstance(starts, list):
            raise ValueError(errmsg)

        if not isinstance(stops, list):
            raise ValueError(errmsg)

        if not isinstance(strands, list):
            raise ValueError(errmsg)

        if not isinstance(scores, list):
            raise ValueError(errmsg)

        if not isinstance(pvalues, list):
            raise ValueError(errmsg)

        if not isinstance(frequencies, list):
            raise ValueError(errmsg)

        if not isinstance(references, list):
            raise ValueError(errmsg)

        self._seqnames = seqnames
        self._seqs = seqs
        self._chroms = chroms
        self._starts = starts
        self._stops = stops
        self._strands = strands
        self._scores = scores
        self._pvalues = pvalues
        self._frequencies = frequencies
        self._references = references


    def get_seqnames(self):
        if not self._seqnames:
            errmsg: str = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._seqnames


    def get_seqs(self):
        if not self._seqs:
            errmsg: str = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._seqs


    def get_chroms(self):
        if not self._chroms:
            errmsg: str = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._chroms


    def get_starts(self):
        if not self._starts:
            errmsg: str = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._starts


    def get_stops(self):
        if not self._stops:
            errmsg: str = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._stops


    def get_strands(self):
        if not self._strands:
            errmsg: str = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._strands


    def get_scores(self):
        if not self._scores:
            errmsg: str = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._scores


    def get_pvalues(self):
        if not self._pvalues:
            errmsg: str = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._pvalues


    def get_frequencies(self):
        if not self._frequencies:
            errmsg: str = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._frequencies


    def get_references(self):
        if not self._references:
            errmsg: str = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._references

# end of ResultTmp

