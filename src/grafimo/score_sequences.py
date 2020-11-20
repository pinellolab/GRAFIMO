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
from grafimo.utils import die, printProgressBar, sigint_handler
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



def compute_results(motif: Motif,
                    sequence_loc: str,
                    args_obj: Optional[Findmotif] = None,
                    testmode: Optional[bool] = False
) -> pd.DataFrame:
    """Score all the sequences extracted from the genome variation graph
    in the regions defined in the input BED file.

    To score the sequences is used the scaled motif scoring matrix, 
    stored in the input Motif instance.

    To each score is assigned a P-value using the P-value matrix, 
    contained in the Motif instance.
    
    Parameters
    ----------
    motif : Motif
        motif data to score sequences
    sequence_loc : str
        path to the intermediate files containing the sequences 
        extracted from the genome variation graph
    args_obj : Findmotif, optional
        container for the arguments needed during the scoring step
    testmode : bool, optional
        flag value manually set used for test purposes

    Returns
    -------
    pandas.DataFrame
        scoring results
    """

    cores:int
    threshold: float
    no_qvalue: bool
    qval_t: bool
    no_reverse: bool
    recomb: bool
    verbose: bool
    errmsg: str

    if not isinstance(sequence_loc, str):
        errmsg = ''.join(["\n\nERROR: unable to locate extracted sequences in ", 
                          sequence_loc])
        raise FileNotFoundError(errmsg)

    if not isinstance(motif, Motif):
        errmsg = "\n\nERROR: the given motif is not an instance of Motif"
        raise ValueError(errmsg)

    if not testmode:
        if not isinstance(args_obj, Findmotif):
            errmsg = "\n\nERROR: unrecognized argument object type"
            raise ValueError(errmsg)

    if not testmode:
        cores = args_obj.get_cores()
        threshold = args_obj.get_threshold()
        no_qvalue = args_obj.get_no_qvalue()
        qval_t = args_obj.get_qvalueT()
        no_reverse = args_obj.get_no_reverse()
        recomb = args_obj.get_recomb()
        verbose = args_obj.get_verbose()
    else:
        cores = 1
        threshold = 1
        recomb = True
        no_qvalue = False
        qval_t = False
        no_reverse = False
        verbose = False

    assert threshold > 0
    assert threshold <= 1
    assert cores >= 1

    print_scoring_msg(no_reverse, motif)

    cwd: str = os.getcwd()
    os.chdir(sequence_loc)

    manager: SyncManager = mp.Manager()
    # results
    return_dict: DictProxy = manager.dict()
    # scanned nucleotides
    scanned_nucs_dict: DictProxy = manager.dict()
    # scanned sequences  
    scanned_seqs_dict: DictProxy = manager.dict()  

    # get all tmp files containing sequences
    sequences: List[str] = glob.glob('*.tsv')  
    if len(sequences) < cores:
        cores = len(sequences)
    # split the sequence set in no. cores chunks
    sequences_split: List[str] = np.array_split(sequences, cores)  

    jobs = list()  # jobs list
    proc_finished: int = 0 

    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    # overwrite the default SIGINT handler to exit gracefully
    # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
    signal.signal(signal.SIGINT, original_sigint_handler)  
    

    if verbose:
        start_s: float = time.time()

    try:

        # compute results in parallel
        for i in range(cores):
            p = mp.Process(
                target=score_seqs, args=(
                    sequences_split[i], motif, no_reverse, return_dict, 
                    scanned_seqs_dict, scanned_nucs_dict, i
                    )
                )
            jobs.append(p)
            p.start()  
        # end for

        # to print 0%, otherwise start from  % as first chunk id already completed completed
        printProgressBar(proc_finished, cores, prefix='Progress:',
                         suffix='Complete', length=50)
        for job in jobs:
            job.join()  # sync point
            proc_finished += 1
            printProgressBar(proc_finished, cores, prefix='Progress:',
                             suffix='Complete', length=50)
        # end for
    
    except KeyboardInterrupt:
        sigint_handler()
        sys.exit(2)

    else:
        if verbose:
            end_s: float = time.time()
            print(
                "Scored all sequences in %.2fs" % (end_s - start_s)
            )

        else:
            pass # all was OK, go to the next instruction
    
    # end try

    os.chdir(cwd) 

    if not testmode:
        cmd: str = "rm -rf {0}".format(sequence_loc)
        code: int = subprocess.call(cmd, shell=True)

        if code != 0:
            errmsg = "\n\nERROR: an error occurred while running %s" % cmd
            raise SubprocessError(errmsg)
    

    if verbose:
        start_df: str = time.time()

    # recover all analysis results and summarize them in a single 
    # data structure
    seqnames: List[str] = list()
    seqs: List[str] = list()
    chroms: List[str] = list()
    starts: List[int] = list()
    stops: List[int] = list()
    strands: List[str] = list()
    scores: List[np.double] = list()
    pvalues: List[np.double] = list()
    frequencies: List[int] = list()
    references: List[str] = list()

    seqs_scanned: int = 0
    nucs_scanned: int = 0

    for key in return_dict.keys():

        assert isinstance(return_dict[key], ResultTmp)

        seqnames += return_dict[key].get_seqnames()
        seqs += return_dict[key].get_seqs()
        chroms += return_dict[key].get_chroms()
        starts += return_dict[key].get_starts()
        stops += return_dict[key].get_stops()
        strands += return_dict[key].get_strands()
        scores += return_dict[key].get_scores()
        pvalues += return_dict[key].get_pvalues()
        frequencies += return_dict[key].get_frequencies()
        references += return_dict[key].get_references()

        # compute the total number of scanned sequences and nucleotides
        seqs_scanned += scanned_seqs_dict[key]  # the keys are the same as return_dict
        nucs_scanned += scanned_nucs_dict[key]  # the keys are the same as return_dict
    # end for

    qvalues: List[np.double]
    # compute the q-values
    if no_qvalue:
        qvalues = list()  # empty list -> not computed
    else:
        qvalues = compute_qvalues(pvalues)
    # end if

    print("Scanned sequences:", seqs_scanned)
    print("Scanned nucleotides:", nucs_scanned)

    # summarize results in a pandas DF
    finaldf: pd.DataFrame = build_df(motif, seqnames, starts, stops, strands, 
                                     scores, pvalues, qvalues, seqs, frequencies, 
                                     references, threshold, qval_t, no_qvalue, 
                                     recomb)

    if verbose:
        end_df: float = time.time()
        print("\nResults summary built in %.2fs" % (end_df - start_df))
    
    return finaldf

# end of compute_results()


def score_seqs(sequences: List[str],
               motif: Motif,
               no_reverse: bool,
               return_dict: DictProxy,
               scanned_seqs_dict: DictProxy,
               scanned_nucs_dict: DictProxy,
               pid: int
) -> None:
    """Score the retrieved sequences using motif scoring matrix data.

    The partial results are stored in a dictionary with current
    process ID as key for the current entry.

    The different entries will be merged at the end, to obtian the 
    final results report.

    Parameters
    ----------
    sequences : list
        sequences to score
    motif : Motif
        motif object containing thescoring matrix and the P-value matrix
    no_reverse:
        if True only the sequences belonging to the forward strand will
        be scored
    return_dict : multiprocessing.managers.DictProxy
        dictionary where the current chunk of results will be
        stored
    scanned_seqs_dict : mp.managers.DictProxy
        dictionary storing the number of sequences scanned in each 
        sequence chunk
    scanned_nucs_dict : mp.managers.DictProxy
        dictionary storing the number of nucleotides scanned in each
        sequence chunk
    """

    try:
        # get Motif attributes to score sequences
        score_matrix: np.ndarray = motif.getMotif_scoreMatrix()
        pval_mat: np.array = motif.getMotif_pval_mat()
        min_score: int = motif.getMin_val()
        scale: int = motif.getScale()
        width: int = motif.getWidth()
        offset: np.double = motif.getOffset()

        # initialize lists where results will be stored
        seqs: List[str] = list()
        scores: List[np.double] = list()
        pvalues: List[np.double] = list()
        seqnames: List[str] = list()
        chroms: List[str] = list()
        starts: List[int] = list()
        stops: List[int] = list()
        strands: List[str] = list()
        frequencies: List[int] = list()
        references: List[str] = list()

        seqs_scanned: int = 0  # counter for scanned sequences 

        width: int = motif.getWidth()

        for s in sequences:
            with open(s, mode='r') as raw_sequences:
                for line in raw_sequences:
                    data = line.split('\t')  

                    strand = data[2][-1]

                    if no_reverse:  # score only the fw strand

                        if strand == '+':

                            # read the final values
                            seq = data[1]
                            seqname = ''.join(['chr', data[0]])
                            chrom = seqname.split(':')[0]
                            start = data[2].split(':')[1]
                            start = start[:-1]
                            stop = data[3].split(':')[1]
                            stop = stop[:-1]
                            freq = data[4]
                            ref = data[5]
                            score, pvalue = compute_score_seq(seq, score_matrix, 
                                                              pval_mat, min_score, 
                                                              scale, width, offset)
                            seqs_scanned += 1
                            seqs.append(seq)
                            scores.append(score)
                            pvalues.append(pvalue)
                            seqnames.append(seqname)
                            chroms.append(chrom)
                            starts.append(start)
                            stops.append(stop)
                            strands.append(strand)
                            frequencies.append(freq)
                            
                            # fix indels reference report bug
                            distance: int = np.abs(int(stop) - int(start))
                            if (ref == "ref" and distance != width):
                                ref = "non.ref"

                            references.append(ref)
                        # end if

                    else:  # score both fw and reverse strands

                        seq = data[1]
                        seqname = ''.join(['chr', data[0]])
                        chrom = seqname.split(':')[0]
                        start = data[2].split(':')[1]
                        start = start[:-1]
                        stop = data[3].split(':')[1]
                        stop = stop[:-1]
                        freq = data[4]
                        ref = data[5]
                        score, pvalue = compute_score_seq(seq, score_matrix, 
                                                          pval_mat, min_score, 
                                                          scale, width, offset)
                        seqs_scanned += 1
                        seqs.append(seq)
                        scores.append(score)
                        pvalues.append(pvalue)
                        seqnames.append(seqname)
                        chroms.append(chrom)
                        starts.append(start)
                        stops.append(stop)
                        strands.append(strand)
                        frequencies.append(freq)

                        # fix indels reference report bug
                        distance: int = np.abs(int(stop) - int(start))
                        if (ref == "ref" and distance != width):
                            ref = "non.ref"

                        references.append(ref)
                    # end if
                # end for
            # end open
        # end for

    except KeyboardInterrupt:
        pass

    else:

        res_tmp = ResultTmp(seqnames, seqs, chroms, starts, stops, strands, 
                            scores, pvalues, frequencies, references)

        return_dict[pid] = res_tmp
        scanned_seqs_dict[pid] = seqs_scanned
        scanned_nucs_dict[pid] = seqs_scanned * width 
    # end try

# end of score_seqs()


@jit(nopython=True)
def compute_score_seq(seq: str,
                      score_matrix: np.ndarray,
                      pval_mat: np.array,
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
    pval_mat : numpy.array
        motif P-value matrix
    min_score : int
        minimum score within the scoring matrix
    scale : int
        scaling factor
    width : int
        motif width
    offset : numpy.double
        scaling offset

    Returns
    -------
    numpy.double
        sequence log-odds score
    numpy.double
        sequence P-value
    """

    score: int = 0

    seq_len: int = len(seq)
    assert seq_len == width

    # score the current sequence
    for i in range(width):
        nuc: str = seq[i]

        if nuc == 'N':
            score = min_score
            break  # we don't go further

        elif nuc == 'A' or nuc == 'a':
            score += score_matrix[0, i]
        elif nuc == 'C' or nuc == 'c':
            score += score_matrix[1, i]
        elif nuc == 'G' or nuc == 'g':
            score += score_matrix[2, i]
        elif nuc == 'T' or nuc == 't':
            score += score_matrix[3, i]
        # end if
    # end for

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


def compute_qvalues(pvalues: List[np.double]) -> List[np.double]:
    """Compute q-values for a given list of P-values.

    The q-values are obtained using the Benjamini-Hochberg method.

    Parameters
    ----------
    pvalues : list
        list of P-values

    Returns
    -------
    list
        list of q-values
    """

    if not isinstance(pvalues, list):
        errmsg: str = "\n\nERROR: P-values must be in a list"
        raise ValueException(errmsg)

    print("\nComputing q-values...\n")

    # use Benjamini-Hochberg procedure to correct P-values
    mt_obj: Tuple[np.ndarray, np.ndarray, np.double, float]
    mt_obj = multipletests(pvalues, method="fdr_bh")
    qvalues: List[float] = list(mt_obj[1])

    return qvalues

# end of compute_qvalues()


def build_df(motif: Motif,
             seqnames: List[str],
             starts: List[int],
             stops: List[int],
             strands: List[str],
             scores: List[np.double],
             pvalues: List[np.double],
             qvalues: List[np.double],
             sequences: List[str],
             frequencies: List[int],
             references: List[str],
             threshold: float,
             qval_t: bool,
             no_qvalue: bool, 
             recomb: bool
) -> pd.DataFrame:
    """Build the results summary report. The results are stored in a 
    pandas DataFrame object.

    The motif occurrence candidates are filtered applying a threshold on
    the P-value or on the q-value.

    The remaining entries are reported in the final results.

    Parameters
    ----------
    motif : Motif
        Motif object
    seqnames : list
        sequence names
    starts : list
        starting coordinates
    stops : list
        stopping coordinates
    strands : list
        DNA strands
    pvalues: list
        P-values
    qvalues : list
        q-values
    sequences : list
        sequences
    references : list
        flag values stating if the sequences contain genomi variants
    threshold : float
        threshold to apply on P-values or q-values in order to filter
        the motif occurrence candidates to report
    qval_t : bool
        if True the threshold will be applied on q-values rather on
        P-values
    no_qvalue:
        if True the q-values have not been computed
    recomb : bool
        if True will be reported also sequences which can be built with 
        the given set of genomic variants but do not appear in the
        available samples haplotypes
    
    Returns
    -------
    pandas.DataFrame
        final results report
    """

    errmsg: str = "\n\nERROR: unknown data-type for motif"
    if not isinstance(motif, Motif):
        raise ValueException(errmsg)

    if not isinstance(seqnames, list):
        raise ValueException(errmsg)

    if not isinstance(starts, list):
        raise ValueException(errmsg)

    if not isinstance(stops, list):
        raise ValueException(errmsg)

    if not isinstance(strands, list):
        raise ValueException(errmsg)

    if not isinstance(pvalues, list):
        raise ValueException(errmsg)

    if not isinstance(qvalues, list):
        raise ValueException(errmsg)

    if not isinstance(sequences, list):
        raise ValueException(errmsg)

    if not isinstance(references, list):
        raise ValueException(errmsg)

    if not isinstance(references, list):
        raise ValueException(errmsg)

    if not isinstance(qval_t, bool):
        raise ValueException(errmsg)

    if not isinstance(no_qvalue, bool):
        raise ValueException(errmsg)

    if not isinstance(recomb, bool):
        raise ValueException(errmsg)

    lst_len: int = len(seqnames)

    assert len(starts) == lst_len
    assert len(stops) == lst_len
    assert len(strands) == lst_len
    assert len(scores) == lst_len
    assert len(pvalues) == lst_len
    assert len(sequences) == lst_len
    assert len(frequencies) == lst_len
    assert len(references) == lst_len

    # check if we want also the q-values
    if not no_qvalue: 
        assert len(qvalues) == lst_len

    # apply the threshold on the q-values rather than on P-values
    if qval_t:  
        assert (not no_qvalue)
        assert len(qvalues) > 0 

    seqnames_thresh: List[str] = list() 
    starts_thresh: List[int] = list()
    ends_thresh: List[int] = list()
    strands_thresh: List[str] = list() 
    scores_thresh: List[np.double] = list()
    pvalues_thresh: List[np.double] = list()
    sequences_thresh: List[str] = list()
    frequencies_thresh: List[int] = list()
    references_thresh: List[str] = list()

    if not no_qvalue:
        qvalues_thresh: List[np.double] = list()

    for i in range(lst_len):

        # ignore binding site candidates which does not appear in any sample
        # if not required by tyhe user to analyze them
        if not recomb and int(frequencies[i]) == 0:
                continue

        if not qval_t:  # apply threshold on P-values
            pvalue: np.double = pvalues[i]

            if pvalue < threshold:
                # only the sequences with a P-value under the threshold survive
                seqnames_thresh.append(seqnames[i])
                starts_thresh.append(starts[i])
                ends_thresh.append(stops[i])
                strands_thresh.append(strands[i])
                scores_thresh.append(scores[i])
                pvalues_thresh.append(pvalues[i])
                sequences_thresh.append(sequences[i])
                frequencies_thresh.append(frequencies[i])
                references_thresh.append(references[i])

                if not no_qvalue:
                    qvalues_thresh.append(qvalues[i])
            # end if

        else:  # apply threshold on q-values
            qvalue: np.double = qvalues[i]

            if qvalue < threshold:

                # only the sequences with a q-value under the threshold survive
                seqnames_thresh.append(seqnames[i])
                starts_thresh.append(starts[i])
                ends_thresh.append(stops[i])
                strands_thresh.append(strands[i])
                scores_thresh.append(scores[i])
                pvalues_thresh.append(pvalues[i])
                sequences_thresh.append(sequences[i])
                frequencies_thresh.append(frequencies[i])
                references_thresh.append(references[i])

                # the last control statement, in the if, in this case is not
                # necessary (we must have the q-values)
                # otherwise we should not be here
                qvalues_thresh.append(qvalues[i])
            # end if
        # end if
    # end for

    df_len: int = len(seqnames_thresh)

    # TF's name and ID list
    motif_ids: List[str] = [motif.getMotifID()] * df_len
    motif_names: List[str] = [motif.getMotifName()] * df_len

    df = pd.DataFrame()
    df['motif_id'] = motif_ids
    df['motif_alt_id'] = motif_names
    df['sequence_name'] = seqnames_thresh
    df['start'] = starts_thresh
    df['stop'] = ends_thresh
    df['strand'] = strands_thresh
    df['score'] = scores_thresh
    df['p-value'] = pvalues_thresh

    # add the q-values to the final data frame if they have been computed
    if not no_qvalue:
        df['q-value'] = qvalues_thresh

    # finish to build the data frame
    df['matched_sequence'] = sequences_thresh
    df['haplotype_frequency'] = frequencies_thresh 
    df['reference'] = references_thresh

    # sort entries by p-value
    df = df.sort_values(['p-value'], ascending=True)

    # reindex the data frame in order to have indexes in range [1, (df_len + 1)]
    df.index = list(range(1, (df_len + 1)))

    return df

# end of build_df()


def print_scoring_msg(no_reverse, motif):
    """Print a message to display on terminal during scoring step of
    GRAFIMO analysis.

    Parameters
    ----------
    no_reverse : bool
        if True will be considered only the forward DNA strand
    motif : Motif
        Motif object
    """

    if not isinstance(motif, Motif):
        errmsg: str = '\n\nERROR: The given motif is not an instance of Motif'
        raise ValueException(errmsg)

    motif_id: str = motif.getMotifID()
    fw_id: str = ''.join(['+', motif_id])

    # we take into account also the reverse complement
    if not no_reverse:
        rev_id: str = ''.join(['-', motif_id])

    print('\nScoring hits for motif', fw_id)

    # if we score also the reverse complement
    if not no_reverse:
        print('Scoring hits for motif', rev_id, end="\n\n")

# end of print_scoring_msg()

