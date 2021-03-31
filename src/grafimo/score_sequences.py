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
from grafimo.GRAFIMOException import WrongPathException, ValueException, \
    SubprocessError
from grafimo.utils import die, printProgressBar, sigint_handler, exception_handler

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

    The potential motif occurrences are scored using the scaled scoring matrix.
    The scaled values are then used to retrieve the corresponding P-value.

    ...
    
    Parameters
    ----------
    motif : Motif
        motif object
    sequence_loc : str
        path to sequences extracted
    debug : bool
        trace the full error stack
    args_obj : Findmotif, optional
        commandline arguments container
    testmode : bool, optional
        test (manually set)

    Returns
    -------
    pandas.DataFrame
        results
    """

    if not isinstance(motif, Motif):
        errmsg = "Expected Motif, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif).__name__), debug)
    if not isinstance(sequence_loc, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(sequence_loc).__name__), debug)
    if not os.path.isdir(sequence_loc):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(sequence_loc), debug)
    if not testmode:
        if not isinstance(args_obj, Findmotif):
            errmsg = "Expected Findmotif, got {}.\n"
            exception_handler(TypeError, errmsg.format(type(args_obj).__name__), debug)

    if not testmode:
        cores: int = args_obj.cores
        threshold: float = args_obj.threshold
        no_qvalue: bool = args_obj.noqvalue
        qval_t: bool = args_obj.qvalueT
        no_reverse: bool = args_obj.noreverse
        recomb: bool = args_obj.recomb
        verbose: bool = args_obj.verbose
    else:  # pytest - during normal execution we should never go here
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
    cwd: str = os.getcwd()
    os.chdir(sequence_loc)
    manager: SyncManager = mp.Manager()
    return_dict: DictProxy = manager.dict()  # results
    scanned_nucs_dict: DictProxy = manager.dict()  # scanned nucleotides 
    scanned_seqs_dict: DictProxy = manager.dict()  # scanned sequences 
    sequences: List[str] = glob.glob('*.tsv')  # sequences
    if len(sequences) < cores: cores = len(sequences)
    # split the sequence set in no. cores chunks
    sequences_split: List[str] = np.array_split(sequences, cores)  
    jobs = list()  # jobs list
    proc_finished: int = 0 
    # overwrite the default SIGINT handler to exit gracefully
    # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)  
    if verbose: start_s: float = time.time()
    try:
        for i in range(cores):
            p = mp.Process(
                target=score_seqs, args=(
                    sequences_split[i], motif, no_reverse, return_dict, 
                    scanned_seqs_dict, scanned_nucs_dict, i, debug
                )
            )
            jobs.append(p)
            p.start()  
        # to print 0%, otherwise start from % as first chunk id already completed completed
        printProgressBar(proc_finished, cores, prefix='Progress:', suffix='Complete', length=50)
        for job in jobs:
            job.join()  # sync point
            proc_finished += 1
            printProgressBar(
                proc_finished, cores, prefix='Progress:', suffix='Complete', length=50
            )
    except KeyboardInterrupt:
        sigint_handler()
        die(2)
    else:
        if verbose:
            end_s: float = time.time()
            print( "Scored all sequences in %.2fs" % (end_s - start_s))
    os.chdir(cwd) 
    if not testmode:
        cmd: str = "rm -rf {}".format(sequence_loc)
        code: int = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = "An error occurred while executing {}.\n"
            exception_handler(SubprocessError, errmsg.format(cmd), debug)
    if verbose: start_df: str = time.time()
    # recover all analysis results and summarize them in a single 
    # data structure
    seqs_scanned: int = 0
    nucs_scanned: int = 0
    summary = ResultTmp()
    for key in return_dict.keys():
        partialres = return_dict[key]
        summary.append_list(
            partialres[0], partialres[1], partialres[2], partialres[3], 
            partialres[4], partialres[5], partialres[6], partialres[7], 
            partialres[8], partialres[9]
        )
        seqs_scanned += scanned_seqs_dict[key]
        nucs_scanned += scanned_nucs_dict[key] 
    if summary.isempty():
        errmsg = "No result retrieved. Unable to proceed. Are you using the correct VGs and searching on the right chromosomes?\n"
        exception_handler(ValueError, errmsg, debug)
    # compute the q-values
    if not no_qvalue:
        if verbose: start_q = time.time()
        qvalues = compute_qvalues(summary.pvalues, debug)
        summary.add_qvalues(qvalues)
        if verbose:
            end_q = time.time()
            print("Q-values computed in %.2fs." % (end_q - start_q))
    print("Scanned sequences:\t{}".format(seqs_scanned))
    print("Scanned nucleotides:\t{}".format(nucs_scanned))
    # summarize results in a pandas DataFrame
    finaldf = summary.to_df(motif, threshold, qval_t, recomb, ignore_qvals=no_qvalue)
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
    sequences : list
        sequences
    motif : Motif
        motif object 
    noreverse : bool
        skip reverse strand
    return_dict : multiprocessing.managers.DictProxy
        result dictionary
    scanned_seqs_dict : mp.managers.DictProxy
        number of scanned sequences dictionary
    scanned_nucs_dict : mp.managers.DictProxy
        number of scanned nucleotides dictionary
    pid : int
        process id
    chroms_prefix : str
        chromosome name prefix
    namemap : dict
        chromosome namemap
    debug : bool
        trace the full error stack
    """

    try:
        assert motif.isScaled
        score_matrix: np.ndarray = motif.scoreMatrix
        pval_mat: np.array = motif.pvalMatrix
        min_score: int = motif.minval
        scale: int = motif.scale
        width: int = motif.width
        offset: np.double = motif.offset
        nucsmap: dict = motif.nucsmap

        #restmp = ResultTmp()
        restmp = [[], [], [], [], [], [], [], [], [], []]
        seqs_scanned: int = 0  # counter for scanned sequences 
        for s in sequences:
            ifstream = open(s, mode="r")
            while True:
                line = ifstream.readline()
                if not line: break  # EOF
                data = line.strip().split()
                strand = data[2][-1]
                if noreverse and strand == "-": continue
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
                        seq, score_matrix, pval_mat, min_score, scale, width, offset
                    )
                    seqs_scanned += 1
                    # fix indel reference report bug
                    distance: int = np.abs(stop - start)
                    if (ref == "ref" and distance != width): ref = "non.ref"
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
    """Assign to a DNA sequence a log-odds score based on motif scoring
    matrix. To each score is assigned a corresponding P-value.

    Parameters
    ----------
    seq : str
        sequence
    score_matrix : numpy.ndarray 
        motif scaled scoring matrix
    pval_mat : numpy.ndarray
        motif P-value matrix
    min_score : int
        minimum score 
    scale : int
        scaling factor
    width : int
        motif width
    offset : numpy.double
        scaling offset

    Returns
    -------
    numpy.double
        log-odds score
    numpy.double
        P-value
    """

    score: int = 0
    for i in range(width):
        nuc: str = seq[i]
        if nuc == 'N':
            score = min_score
            break  # we don't go further
        if nuc.upper() == "A": nucidx = 0
        elif nuc.upper() == "C": nucidx = 1
        elif nuc.upper() == "G": nucidx = 2
        elif nuc.upper() == "T": nucidx = 3
        score += score_matrix[nucidx, i]
    assert score >= min_score
    # get the p-value for the obtained score
    tot = pval_mat.sum()
    pvalue: np.double = (pval_mat[score:].sum()) / tot
    # retrieve the log-likelihood score
    logodds: np.double = (score / scale) + (width * offset)
    score = logodds
    assert (pvalue > 0 and pvalue <= 1)
    return score, pvalue

# end of compute_score_seq()


def compute_qvalues(pvalues: List[np.double], debug: bool) -> List[np.double]:
    """Corrects P-values with False Discovery Rate Benjamini-Hochberg procedure.

    ...

    Parameters
    ----------
    pvalues : list
        P-values
    debug : bool
        trace the full error stack

    Returns
    -------
    list
        corrected P-values (q-values)
    """

    if not isinstance(pvalues, list):
        errmsg = "Expected list, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(pvalues).__name__), debug)

    print("\nComputing q-values...\n")
    # use Benjamini-Hochberg procedure to correct P-values
    mt_obj = multipletests(pvalues, method="fdr_bh")
    qvalues: List[float] = list(mt_obj[1])

    return qvalues

# end of compute_qvalues()


def print_scoring_msg(motif: Motif, noreverse: bool, debug: bool):
    """Message printed when scoring procedure begins.

    ...

    Parameters
    ----------
    motif : Motif
        motif object
    noreverse : bool
        skip reverse strand sequences
    debug : bool
        trace the full error stack
    """

    if not isinstance(motif, Motif):
        errmsg = "Expected Motif, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif).__name__), debug)
    if not isinstance(noreverse, bool):
        errmsg = "Expected bool, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(noreverse).__name__), debug)

    fw_id: str = "".join(["+", motif.motifID])
    if not noreverse: rev_id: str = "".join(["-", motif.motifID])
    msg = "Scoring hits for motif {}."
    print(msg.format(fw_id))
    if not noreverse:
        print(msg.format(rev_id), end="\n\n")

# end of print_scoring_msg()

