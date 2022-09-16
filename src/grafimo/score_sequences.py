"""Functions to score the motif occurrence candidates found scanning
the genome variation graph in the given genomic regions.

To score the motif occurrences is used the scaled scoring matrix 
previously computed from motif PWM data. The resulting scaled scores are
then converted in the corresponding log-odds values. To each score
is given a P-value and if requested a q-value.

The occurrences retrieved in the results are only those surviving the 
threshold put by the user on the P-value or on the q-value (1e-4 by
default on both)

The running time of the algorithm used is bounded by O(n^2).
"""


from grafimo.motif import Motif
from grafimo.workflow import Findmotif
from grafimo.resultsTmp import ResultTmp
from grafimo.utils import (
    die, 
    print_progress_bar, 
    sigint_handler, 
    exception_handler
)

from typing import List, Optional, Dict, Tuple
from statsmodels.stats.multitest import multipletests
from multiprocessing.managers import DictProxy, SyncManager
from numba import jit

import multiprocessing as mp
import pandas as pd
import numpy as np

import subprocess
import signal
import time
import glob
import sys
import os


def compute_results(
    motif: Motif,
    sequence_loc: str,
    debug: bool,
    args_obj: Optional[Findmotif] = None,
    testmode: Optional[bool] = False,
) -> pd.DataFrame:
    """Score the sequences extracted from the genome variation graph.

    The potential motif occurrences are scored using the scaled scoring 
    matrix. The scaled values are then used to retrieve the corresponding 
    P-value.

    ...
    
    Parameters
    ----------
    motif : Motif
        Motif 
    sequence_loc : str
        Path to sequences location
    debug : bool
        Trace the full error stack
    args_obj : Findmotif, optional
        Commandline arguments container
    testmode : bool, optional
        Test mode (manually set)

    Returns
    -------
    pandas.DataFrame
        Results table
    """

    if not isinstance(motif, Motif):
        errmsg = f"Expected {type(Motif).__name__}, got {type(motif).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(sequence_loc, str):
        errmsg = f"Expected {str.__name__}, got {type(sequence_loc).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not os.path.isdir(sequence_loc):
        errmsg = f"Unable to locate {sequence_loc}.\n"
        exception_handler(FileNotFoundError, errmsg, debug)
    if not testmode:
        if not isinstance(args_obj, Findmotif):
            errmsg = f"Expected {type(Findmotif).__name__}, got {type(args_obj).__name__}.\n"
            exception_handler(TypeError, errmsg, debug)
    # read input arguments
    if not testmode:
        cores = args_obj.cores
        threshold = args_obj.threshold
        no_qvalue = args_obj.noqvalue
        qval_t = args_obj.qvalueT
        no_reverse = args_obj.noreverse
        recomb = args_obj.recomb
        verbose = args_obj.verbose
    else:  # for pytest ONLY - normally we should NEVER go here
        cores = 1
        threshold = float(1)
        recomb = True
        no_qvalue = False
        qval_t = False
        no_reverse = False
        verbose = False
    assert threshold > 0 and threshold <= 1
    assert cores >= 1
    print_scoring_msg(motif, no_reverse, debug)
    cwd = os.getcwd()
    # recover the correct potential motif matches
    sequence_loc = os.path.join(sequence_loc, f"width_{motif.width}")
    os.chdir(sequence_loc)
    manager = mp.Manager()  # manager to recover multiprocessing results
    return_dict = manager.dict()  # results
    scanned_nucs_dict = manager.dict()  # scanned nucleotides 
    scanned_seqs_dict = manager.dict()  # scanned sequences 
    sequences = glob.glob("*.tsv")  # sequences
    if len(sequences) < cores: 
        cores = len(sequences)  # avoid wasting cores 
    # split the sequence set in no. cores chunks
    sequences_split = np.array_split(sequences, cores)  
    jobs = []  # jobs list
    proc_finished = 0 
    # overwrite the default SIGINT handler to exit gracefully
    # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)  
    if verbose: 
        start_s = time.time()
    try:
        for i in range(cores):
            p = mp.Process(
                target=score_seqs, args=(
                    sequences_split[i], 
                    motif, 
                    no_reverse, 
                    return_dict, 
                    scanned_seqs_dict, 
                    scanned_nucs_dict, 
                    i, 
                    debug
                )
            )
            jobs.append(p)
            p.start()  
        # to print 0%, otherwise start from % as first chunk id already completed completed
        print_progress_bar(
            proc_finished, cores, prefix="Progress:", suffix="Complete", length=50
        )
        for job in jobs:
            job.join()  # sync point
            proc_finished += 1
            print_progress_bar(
                proc_finished, cores, prefix="Progress:", suffix="Complete", length=50
            )
    except KeyboardInterrupt:
        sigint_handler()
        die(2)
    else:
        if verbose:
            end_s = time.time()
            print("Sequences scored in %.2fs" % (end_s - start_s))
    os.chdir(cwd) # go back to source location
    # recover all analysis results and summarize them in a single report
    if verbose: 
        start_df = time.time()
    seqs_scanned = 0
    nucs_scanned = 0
    summary = ResultTmp()
    for key in return_dict.keys():
        partialres = return_dict[key]
        # TODO: improve partial results recovery
        summary.append_list(
            partialres[0], 
            partialres[1], 
            partialres[2], 
            partialres[3], 
            partialres[4], 
            partialres[5], 
            partialres[6], 
            partialres[7], 
            partialres[8], 
            partialres[9]
        )
        seqs_scanned += scanned_seqs_dict[key]
        nucs_scanned += scanned_nucs_dict[key] 
    if summary.isempty():  # no results?
        errmsg = "No result retrieved. Unable to proceed.\n" 
        errmsg += "\nAre you using the correct VGs and searching on the right chromosomes?\n"
        exception_handler(ValueError, errmsg, debug)
    # compute q-values
    if not no_qvalue:
        if verbose: 
            start_q = time.time()
        qvalues = compute_qvalues(summary.pvalues, debug)
        summary.add_qvalues(qvalues)
        if verbose:
            end_q = time.time()
            print("Q-values computed in %.2fs." % (end_q - start_q))
    print(f"Scanned sequences:\t{seqs_scanned}")
    print(f"Scanned nucleotides:\t{nucs_scanned}")
    # summarize results in a pandas DataFrame
    finaldf = summary.to_df(
        motif, threshold, qval_t, recomb, ignore_qvals=no_qvalue
    )
    if verbose:
        end_df: float = time.time()
        print("\nResults summary built in %.2fs" % (end_df - start_df))
    return finaldf

# end of compute_results()


def score_seqs(
    sequences: List[str],
    motif: Motif,
    noreverse: bool,
    return_dict: DictProxy,
    scanned_seqs_dict: DictProxy,
    scanned_nucs_dict: DictProxy,
    pid: int,
    debug: bool
) -> None:
    """Score sequences extracted from genome variation graph(s).
    
    The partial results are stored in a dictionary. The process ids are used as
    keys of the dictionary.

    ...

    Parameters
    ----------
    sequences : List[str]
        Potential motif occurrence sequences
    motif : Motif
        Motif 
    noreverse : bool
        Skip reverse strand
    return_dict : multiprocessing.managers.DictProxy
        Result dictionary
    scanned_seqs_dict : mp.managers.DictProxy
        Number of scanned sequences dictionary
    scanned_nucs_dict : mp.managers.DictProxy
        Number of scanned nucleotides dictionary
    pid : int
        Process id
    chroms_prefix : str
        Chromosome name prefix
    namemap : dict
        Chromosome namemap
    debug : bool
        Trace the full error stack

    Returns
    -------
    None
    """

    try:
        assert motif.is_scaled
        score_matrix = motif.score_matrix
        pval_mat = motif.pval_matrix
        min_score = motif.min_val
        scale = motif.scale
        width = motif.width
        offset = motif.offset

        #restmp = ResultTmp()
        restmp = [[], [], [], [], [], [], [], [], [], []]
        seqs_scanned = 0  # counter for scanned sequences 
        for s in sequences:
            handle = open(s, mode="r")
            while True:
                line = handle.readline()
                if not line: 
                    break  # EOF
                data = line.strip().split()
                strand = data[2][-1]
                if noreverse and strand == "-": 
                    continue
                else:
                    # parse data
                    seqname = data[0]
                    seq = data[1]
                    chrom = seqname.split(":")[0]
                    start = data[2].split(":")[1]
                    start = int(start[:-1])  # remove strand
                    stop = data[3].split(":")[1]
                    stop = int(stop[:-1])  # remove strand
                    freq = data[4]
                    ref = data[5]
                    score, pvalue = compute_score_seq(
                        seq, 
                        score_matrix, 
                        pval_mat, 
                        min_score, 
                        scale, 
                        width, 
                        offset
                    )
                    seqs_scanned += 1
                    # fix indel reference report bug
                    distance = np.abs(stop - start)
                    if (ref == "ref" and distance != width): 
                        ref = "non.ref"
                    #restmp.append(
                    #    seqname, seq, chrom, start, stop, strand, score, pvalue, 
                    #    int(freq), ref
                    #)
                    restmp[0].append(seqname)
                    restmp[1].append(seq)
                    restmp[2].append(chrom)
                    restmp[3].append(start)
                    restmp[4].append(stop)
                    restmp[5].append(strand)
                    restmp[6].append(score)
                    restmp[7].append(pvalue)
                    restmp[8].append(int(freq))
                    restmp[9].append(ref)
    except KeyboardInterrupt:
        pass  # handled above
    return_dict[pid] = restmp
    scanned_seqs_dict[pid] = seqs_scanned
    scanned_nucs_dict[pid] = seqs_scanned * width 

# end of score_seqs()


@jit(nopython=True)
def compute_score_seq(
    seq: str,
    score_matrix: np.ndarray,
    pval_mat: np.ndarray,
    min_score: int,
    scale: int,
    width: int,
    offset: np.double
) -> Tuple[np.double, np.double]:
    """Compute the likelihood score for a DNA sequece to be a potential
    occurrence of the input motif. To each score is assigned the 
    corresponding statistical significance (P-value). The P-values are 
    then corrected in q-values.
    
    ...

    Parameters
    ----------
    seq : str
        DNA sequence
    score_matrix : numpy.ndarray 
        Motif scaled scoring matrix
    pval_mat : numpy.ndarray
        Motif P-value matrix
    min_score : int
        Minimum score in the scaled motif matrix 
    scale : int
        Scaling factor
    width : int
        Motif width
    offset : numpy.double
        Scaling offset

    Returns
    -------
    numpy.double
        Sequence log-odds score
    numpy.double
        Log-odds score P-value
    """

    score = 0
    for i in range(width):
        nuc = seq[i]
        if nuc == "N":
            score = min_score
            break  # we don't go further
        if nuc.upper() == "A": 
            nucidx = 0
        elif nuc.upper() == "C": 
            nucidx = 1
        elif nuc.upper() == "G": 
            nucidx = 2
        elif nuc.upper() == "T": 
            nucidx = 3
        score += score_matrix[nucidx, i]
    assert score >= min_score
    # get the p-value for the obtained score
    tot = pval_mat.sum()
    pvalue = (pval_mat[score:].sum()) / tot
    # retrieve the log-likelihood score
    logodds = (score / scale) + (width * offset)
    score = logodds
    assert (pvalue > 0 and pvalue <= 1)
    return score, pvalue

# end of compute_score_seq()


def compute_qvalues(pvalues: List[np.double], debug: bool) -> List[np.double]:
    """Correct P-values by False Discovery Rate (Benjamini-Hochberg 
    procedure).

    ...

    Parameters
    ----------
    pvalues : List[float]
        P-values
    debug : bool
        Trace the full error stack

    Returns
    -------
    List[float]
        Corrected P-values (q-values)
    """

    if not isinstance(pvalues, list):
        errmsg = f"Expected {list.__name__}, got {type(pvalues).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    print("\nComputing q-values...\n")
    # use Benjamini-Hochberg to correct P-values
    res_fdr = multipletests(pvalues, method="fdr_bh")
    qvalues = list(res_fdr[1])
    assert len(pvalues) == len(qvalues)
    return qvalues

# end of compute_qvalues()


def print_scoring_msg(motif: Motif, noreverse: bool, debug: bool) -> None:
    """Print this message when the scoring step begins.

    ...

    Parameters
    ----------
    motif : Motif
        Motif 
    noreverse : bool
        Skip reverse strand sequences
    debug : bool
        Trace the full error stack

    Returns
    -------
    None
    """

    if not isinstance(motif, Motif):
        errmsg = f"Expected {type(Motif).__name__}, got {type(motif).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(noreverse, bool):
        errmsg = f"Expected {bool.__name__}, got {type(noreverse).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    fw_id = "".join(["+", motif.motif_id])
    if not noreverse: 
        rev_id = "".join(["-", motif.motif_id])
    msg = "Scoring hits for motif {}."
    print(msg.format(fw_id))
    if not noreverse:
        print(msg.format(rev_id), end="\n\n")

# end of print_scoring_msg()

