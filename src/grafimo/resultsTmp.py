"""Definition of ResultTmp class.

The class is used to store intermediate results during the motif 
occurrences candidates scoring step of GRAFIMO analysis.
"""

from grafimo.motif import Motif

from typing import List, Optional

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

    # these errors should never appear --> no need for error formatting
    # can assume that debug mode == True
    def __init__(self):
        self._seqnames = list()
        self._seqs = list()
        self._chroms = list()
        self._starts = list()
        self._stops = list()
        self._strands = list()
        self._scores = list()
        self._pvalues = list()
        self._qvalues = list()
        self._frequencies = list()
        self._references = list()


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


    def append(
        self, 
        seqname: str, 
        seq: str, 
        chrom: str, 
        start: int, 
        stop: int, 
        strand: str, 
        score: float, 
        pvalue: float, 
        freq: int, 
        ref: str
    ) -> None:
        if not isinstance(seqname, str):
            errmsg = f"\n\nERROR: Expected {str.__name__}, got {type(seqname).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(seq, str):
            errmsg = f"\n\nERROR: Expected {str.__name__}, got {type(seq).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(chrom, str):
            errmsg = f"\n\nERROR: Expected {str.__name__}, got {type(chrom).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(start, int):
            errmsg = f"\n\nERROR: Expected {int.__name__}, got {type(start).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(stop, int):
            errmsg = f"\n\nERROR: Expected {int.__name__}, got {type(stop).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(strand, str):
            errmsg = f"\n\nERROR: Expected {str.__name__}, got {type(strand).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(score, float):
            errmsg = f"\n\nERROR: Expected {float.__name__}, got {type(score).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(pvalue, float):
            errmsg = f"\n\nERROR: Expected {float.__name__}, got {type(pvalue).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(freq, int):
            errmsg = f"\n\nERROR: Expected {int.__name__}, got {type(freq).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(ref, str):
            errmsg = f"\n\nERROR: Expected {str.__name__}, got {type(ref).__name__}.\n"
            raise TypeError(errmsg)
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
        seqnames: List[str], 
        seqs: List[str], 
        chroms: List[str], 
        starts: List[int], 
        stops: List[int], 
        strands: List[str], 
        scores: List[float], 
        pvalues: List[float], 
        frequencies: List[int], 
        references: List[str]
    ) -> None:
        if not isinstance(seqnames, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(seqnames).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(seqs, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(seqs).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(chroms, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(chroms).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(starts, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(starts).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(stops, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(stops).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(strands, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(strands).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(scores, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(scores).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(pvalues, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(pvalues).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(frequencies, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(frequencies).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(references, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(references).__name__}.\n"
            raise TypeError(errmsg) 
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

    
    def add_qvalues(self, qvalues: List[float]) -> None:
        if not isinstance(qvalues, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(qvalues).__name__}.\n"
            raise TypeError(errmsg)
        self._qvalues = qvalues

    
    def to_df(
        self, 
        motif: Motif, 
        threshold: float, 
        qvalt: bool, 
        recomb: bool, 
        ignore_qvals: Optional[bool] = False
    ) -> pd.DataFrame:
        if not isinstance(motif, Motif):
            errmsg = f"\n\nERROR: Expected {type(Motif).__name__}, got {type(motif).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(threshold, float):
            errmsg = f"\n\nERROR: Expected {float.__name__}, got {type(threshold).__name__}.\n"
            raise TypeError(errmsg)
        if threshold <= 0 or threshold > 1:
            errmsg = "\n\nERROR: The threshold must be between 0 and 1.\n"
            raise ValueError(errmsg)
        if not isinstance(qvalt, bool):
            errmsg = f"\n\nERROR: Expected {bool.__name__}, got {type(qvalt).__name__}.\n"
            raise ValueError(errmsg)
        if not isinstance(recomb, bool):
            errmsg = f"Expected {bool.__name__}, got {type(recomb).__name__}.\n"
            raise TypeError(errmsg)
        if not isinstance(ignore_qvals, bool):
            errmsg = f"\n\nERROR: Expected {bool.__name__}, got {type(ignore_qvals).__name__}.\n"
            raise TypeError(errmsg)
        if qvalt: 
            assert bool(self._qvalues) and not ignore_qvals
        if ignore_qvals:
            df = pd.DataFrame(
                {
                    "motif_id":[motif.motif_id] * len(self._seqnames),
                    "motif_alt_id":[motif.motif_name] * len(self._seqnames),
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
                    "motif_id":[motif.motif_id] * len(self._seqnames),
                    "motif_alt_id":[motif.motif_name] * len(self._seqnames),
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
        df_thresh = df_thresh.sort_values(["p-value"], ascending=True)
        df_thresh.reset_index(drop=True, inplace=True)
        return df_thresh

    
    def isempty(self) -> bool:
        # q-values can be empty --> ignore them
        # if just one of the mandatory fields is empty we cannot proceed 
        if (
            not self._seqnames or not self._seqs or not self._chroms or 
            not self._starts or not self._stops or not self._strands or 
            not self._scores or not self._pvalues or not self._frequencies or
            not self._references
        ):
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

