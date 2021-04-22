"""Definition of ResultTmp class.

The class is used to store intermediate results during the motif 
occurrences candidates scoring step of GRAFIMO analysis.
"""

from grafimo.motif import Motif

from typing import List

import pandas as pd
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
    # ResultTmp Attributes
    #
    _seqnames: List[str] = list()
    _seqs: List[str] = list()
    _chroms: List[str] = list()
    _starts: List[int] = list()
    _stops: List[int] = list()
    _strands: List[str] = list()
    _scores: List[np.double] = list()
    _pvalues: List[np.double] = list()
    _qvalues: List[np.double] = list()
    _frequencies: List[int] = list()
    _references: List[str] = list()

    #-------------------------------------------------------------------
    # ResultTmp methods
    #

    # these errors should never appear --> no need for error formatting
    # can assume that debug mode == True
    def __init__(self):
        pass


    def __str__(self):
        return "\n".join([
            str(self._seqnames),
            str(self._seqnames),
            str(self._seqs),
            str(self._chroms),
            str(self._starts),
            str(self._stops),
            str(self._strands),
            str(self._scores),
            str(self._pvalues),
            str(self._qvalues),
            str(self._frequencies),
            str(self._references),
        ])


    def append(self, seqname, seq, chrom, start, stop, strand, score, pvalue, freq, ref):
        if not isinstance(seqname, str):
            errmsg = "\n\nERROR: Expected str, got {}.\n"
            raise TypeError(errmsg.format(type(seqname).__name__))
        if not isinstance(seq, str):
            errmsg = "\n\nERROR: Expected str, got {}.\n"
            raise TypeError(errmsg.format(type(seq).__name__))
        if not isinstance(chrom, str):
            errmsg = "\n\nERROR: Expected str, got {}.\n"
            raise TypeError(errmsg.format(type(chrom).__name__))
        if not isinstance(start, int):
            errmsg = "\n\nERROR: Expected int, got {}.\n"
            raise TypeError(errmsg.format(type(start).__name__))
        if not isinstance(stop, int):
            errmsg = "\n\nERROR: Expected int, got {}.\n"
            raise TypeError(errmsg.format(type(stop).__name__))
        if not isinstance(strand, str):
            errmsg = "\n\nERROR: Expected str, got {}.\n"
            raise TypeError(errmsg.format(type(strand).__name__))
        if not isinstance(score, float):
            errmsg = "\n\nERROR: Expected float, got {}.\n"
            raise TypeError(errmsg.format(type(score).__name__))
        if not isinstance(pvalue, float):
            errmsg = "\n\nERROR: Expected float, got {}.\n"
            raise TypeError(errmsg.format(type(pvalue).__name__))
        if not isinstance(freq, int):
            errmsg = "\n\nERROR: Expected int, got {}.\n"
            raise TypeError(errmsg.format(type(freq).__name__))
        if not isinstance(ref, str):
            errmsg = "\n\nERROR: Expected str, got {}.\n"
            raise TypeError(errmsg.format(type(ref).__name__))

        self._seqnames.append(seqname)
        self._seqs.append(seq)
        self._chroms.append(chrom)
        self._starts.append(start)
        self._stops.append(stop)
        self._strands.append(strand)
        self._scores.append(score)
        self._pvalues.append(pvalue)
        self._frequencies.append(freq)
        self._references.append(ref)


    def size(self):
        assert len(self._seqnames) == len(self._seqs)
        assert len(self._seqnames) == len(self._chroms)
        assert len(self._seqnames) == len(self._starts)
        assert len(self._seqnames) == len(self._stops)
        assert len(self._seqnames) == len(self._strands)
        assert len(self._seqnames) == len(self._scores)
        assert len(self._seqnames) == len(self._pvalues)
        assert len(self._seqnames) == len(self._frequencies)
        assert len(self._seqnames) == len(self._references)

        return len(self._seqnames)

    
    def append_list(
        self, 
        seqnames, 
        seqs, 
        chroms, 
        starts, 
        stops, 
        strands, 
        scores, 
        pvalues, 
        frequencies, 
        references
    ):
        if not isinstance(seqnames, list):
            errmsg = "Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(seqnames).__name__))
        if not isinstance(seqs, list):
            errmsg = "Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(seqs).__name__))
        if not isinstance(chroms, list):
            errmsg = "Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(chroms).__name__))
        if not isinstance(starts, list):
            errmsg = "Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(starts).__name__))
        if not isinstance(stops, list):
            errmsg = "Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(stops).__name__))
        if not isinstance(strands, list):
            errmsg = "Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(strands).__name__))
        if not isinstance(scores, list):
            errmsg = "Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(scores).__name__))
        if not isinstance(pvalues, list):
            errmsg = "Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(pvalues).__name__))
        if not isinstance(frequencies, list):
            errmsg = "Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(frequencies).__name__))
        if not isinstance(references, list):
            errmsg = "Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(references).__name__))
        
        self._seqnames += seqnames
        self._seqs += seqs
        self._chroms += chroms
        self._starts += starts
        self._stops += stops
        self._strands += strands
        self._scores += scores
        self._pvalues += pvalues
        self._frequencies += frequencies
        self._references += references

    
    def add_qvalues(self, qvalues):
        if not isinstance(qvalues, list):
            errmsg = "\n\nERROR: Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(qvalues).__name__))
        self._qvalues = qvalues

    
    def to_df(self, motif, threshold, qvalt, recomb, ignore_qvals=False):
        if not isinstance(motif, Motif):
            errmsg = "\n\nERROR: Expected Motif, got {}.\n"
            raise TypeError(errmsg.format(type(motif).__name__))
        if not isinstance(threshold, float):
            errmsg = "\n\nERROR: Expected float, got {}.\n"
            raise TypeError(errmsg.format(type(threshold).__name__))
        if threshold <= 0 or threshold > 1:
            errmsg = "\n\nERROR: The threshold must be between 0 and 1.\n"
            raise ValueError(errmsg)
        if not isinstance(qvalt, bool):
            errmsg = "\n\nERROR: Expected bool, got {}.\n"
            raise ValueError(errmsg.format(type(qvalt).__name__))
        if not isinstance(recomb, bool):
            errmsg = "Expected bool, got {}.\n"
            raise TypeError(errmsg.format(type(recomb).__name__))
        if not isinstance(ignore_qvals, bool):
            errmsg = "\n\nERROR: Expected bool, got {}.\n"
            raise TypeError(errmsg.format(type(ignore_qvals).__name__))
        
        if qvalt: assert bool(self._qvalues) and not ignore_qvals
        if ignore_qvals:
            df = pd.DataFrame(
                {
                    "motif_id":[motif.motifID] * len(self._seqnames),
                    "motif_alt_id":[motif.motifName] * len(self._seqnames),
                    "sequence_name":self.seqnames,
                    "start":self.starts,
                    "stop":self.stops,
                    "strand":self.strands,
                    "score":self.scores,
                    "p-value":self.pvalues,
                    "matched_sequence":self.seqs,
                    "haplotype_frequency":self.frequencies,
                    "reference":self.references
                }
            )
        else:  # ignore_qvals == False
            df = pd.DataFrame(
                {
                    "motif_id":[motif.motifID] * len(self._seqnames),
                    "motif_alt_id":[motif.motifName] * len(self._seqnames),
                    "sequence_name":self.seqnames,
                    "start":self.starts,
                    "stop":self.stops,
                    "strand":self.strands,
                    "score":self.scores,
                    "p-value":self.pvalues,
                    "q-value":self.qvalues,
                    "matched_sequence":self.seqs,
                    "haplotype_frequency":self.frequencies,
                    "reference":self.references
                }
            )
        # apply threshold
        if qvalt:
            assert bool(self.qvalues)
            df_thresh = df[df["q-value"] < threshold]
        else:
            df_thresh = df[df["p-value"] < threshold]
        # remove not observed recombinants
        if not recomb:
            df_thresh = df_thresh[df_thresh["haplotype_frequency"] > 0]
        # sort by p-value
        df_thresh = df_thresh.sort_values(["p-value"], ascending=True).reset_index(drop=True)

        return df_thresh

    
    def isempty(self):
        # q-values can be empty --> ignore them
        # if just one of the mandatory fields is empty we cannot proceed 
        if (not self._seqnames or not self._seqs or not self._chroms or 
            not self._starts or not self._stops or not self._strands or 
            not self._scores or not self._pvalues or not self._frequencies or
            not self._references):
            return True
        return False


    def _get_seqnames(self):
        if not self._seqnames:
            errmsg = "\n\nERROR: \"self._seqnames\" is empty.\n"
            raise AttributeError(errmsg)
        return self._seqnames
    
    @property
    def seqnames(self):
        return self._get_seqnames()


    def _get_seqs(self):
        if not self._seqs:
            errmsg = "\n\nERROR: \"self._seqs\" is empty.\n"
            raise AttributeError(errmsg)
        return self._seqs

    @property
    def seqs(self):
        return self._get_seqs()


    def _get_chroms(self):
        if not self._chroms:
            errmsg = "\n\nERROR: \"self._chroms\" is empty.\n" 
            raise AttributeError(errmsg)
        return self._chroms

    @property
    def chroms(self):
        return self._get_chroms()


    def _get_starts(self):
        if not self._starts:
            errmsg = "\n\nERROR: \"self._starts\" is empty.\n"
            raise AttributeError(errmsg)
        return self._starts
    
    @property
    def starts(self):
        return self._get_starts()


    def _get_stops(self):
        if not self._stops:
            errmsg = "\n\nERROR: \"self._stps\" is empty.\n"
            raise AttributeError(errmsg)
        return self._stops

    @property
    def stops(self):
        return self._get_stops()


    def _get_strands(self):
        if not self._strands:
            errmsg = "\n\nERROR: \"self._strands\" is empty.\n"
            raise AttributeError(errmsg)
        return self._strands

    @property
    def strands(self):
        return self._get_strands()


    def _get_scores(self):
        if not self._scores:
            errmsg = "\n\nERROR: \"self._scores\" is empty.\n"
            raise AttributeError(errmsg)
        return self._scores

    @property
    def scores(self):
        return self._get_scores()


    def _get_pvalues(self):
        if not self._pvalues:
            errmsg = "\n\nERROR: \"self._pvalues\" is empty.\n"
            raise AttributeError(errmsg)
        return self._pvalues

    @property
    def pvalues(self):
        return self._get_pvalues()


    def _get_qvalues(self):
        if not self._qvalues:
            errmsg = "\n\nERROR: \"self._qvalues\" is empty.\n"
            raise AttributeError(errmsg)
        return self._qvalues

    @property
    def qvalues(self):
        return self._get_qvalues()


    def _get_frequencies(self):
        if not self._frequencies:
            errmsg = "\n\nERROR: \"self._frequencies\" is empty.\n"
            raise AttributeError(errmsg)
        return self._frequencies

    @property
    def frequencies(self):
        return self._get_frequencies()


    def _get_references(self):
        if not self._references:
            errmsg = "\n\nERROR: \"self._references\" is empty.\n"
            raise AttributeError(errmsg)
        return self._references

    @property
    def references(self):
        return self._get_references()

# end of ResultTmp

