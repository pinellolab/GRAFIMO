"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Scoring of the retrieved sequences.

Each sequence o length L, where L is the length of the motif,
retrieved from the regions defined in the BED file, is scored
using the scoring matrix built from the input motif file.

The running time of the algorithm used is bounded by O(n^2).

"""


from grafimo.motif import Motif
from grafimo.workflow import Findmotif
from grafimo.GRAFIMOException import WrongPathException, ValueException, SubprocessError
from grafimo.utils import die, printProgressBar, sigint_handler
from statsmodels.stats.multitest import multipletests
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


class ResultTmp(object):
    """
        Class to store intermediate results of the
        sequence scoring step
    """

    _seqnames = None
    _seqs = None
    _chroms = None
    _starts = None
    _stops = None
    _strands = None
    _scores = None
    _pvalues = None
    _references = None

    def __init__(self,
                 seqnames,
                 seqs,
                 chroms,
                 starts,
                 stops,
                 strands,
                 scores,
                 pvalues,
                 references):

        assert len(seqnames) == len(seqs)
        assert len(seqnames) == len(chroms)
        assert len(seqnames) == len(starts)
        assert len(seqnames) == len(stops)
        assert len(seqnames) == len(strands)
        assert len(seqnames) == len(scores)
        assert len(seqnames) == len(pvalues)
        assert len(seqnames) == len(references)

        if not isinstance(seqnames, list):
            errmsg = "\n\nERROR: unable to store temporary results. Wrong data-type given"
            raise ValueError(errmsg)

        if not isinstance(seqs, list):
            errmsg = "\n\nERROR: unable to store temporary results. Wrong data-type given"
            raise ValueError(errmsg)

        if not isinstance(chroms, list):
            errmsg = "\n\nERROR: unable to store temporary results. Wrong data-type given"
            raise ValueError(errmsg)

        if not isinstance(starts, list):
            errmsg = "\n\nERROR: unable to store temporary results. Wrong data-type given"
            raise ValueError(errmsg)

        if not isinstance(stops, list):
            errmsg = "\n\nERROR: unable to store temporary results. Wrong data-type given"
            raise ValueError(errmsg)

        if not isinstance(strands, list):
            errmsg = "\n\nERROR: unable to store temporary results. Wrong data-type given"
            raise ValueError(errmsg)

        if not isinstance(scores, list):
            errmsg = "\n\nERROR: unable to store temporary results. Wrong data-type given"
            raise ValueError(errmsg)

        if not isinstance(pvalues, list):
            errmsg = "\n\nERROR: unable to store temporary results. Wrong data-type given"
            raise ValueError(errmsg)

        if not isinstance(references, list):
            errmsg = "\n\nERROR: unable to store temporary results. Wrong data-type given"
            raise ValueError(errmsg)

        self._seqnames = seqnames
        self._seqs = seqs
        self._chroms = chroms
        self._starts = starts
        self._stops = stops
        self._strands = strands
        self._scores = scores
        self._pvalues = pvalues
        self._references = references

    def get_seqnames(self):
        if not self._seqnames:
            errmsg = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._seqnames

    def get_seqs(self):
        if not self._seqs:
            errmsg = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._seqs

    def get_chroms(self):
        if not self._chroms:
            errmsg = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._chroms

    def get_starts(self):
        if not self._starts:
            errmsg = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._starts

    def get_stops(self):
        if not self._stops:
            errmsg = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._stops

    def get_strands(self):
        if not self._strands:
            errmsg = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._strands

    def get_scores(self):
        if not self._scores:
            errmsg = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._scores

    def get_pvalues(self):
        if not self._pvalues:
            errmsg = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._pvalues

    def get_references(self):
        if not self._references:
            errmsg = "\n\nERROR: attempting to access an empty attribute"
            raise ValueException(errmsg)

        return self._references

# end of ResultTmp


def compute_results(motif,
                    sequence_loc,
                    args_obj):
    """
        Score all the sequences extracted from regions defined in the
        input BED file.
        To score sequences is used the processed input motif.
        The results are then stored in a pandas DataFrame
        ----
        Parameters:
            motif (Motif) : processed motif, used to score sequences
            sequence_loc (str) : path to temporary files storing sequences extracted
                                    during the previous step
            args_obj (Findmotif) : arguments used during the sequnece scoring step
        ----
        Returns:
            finaldf (pd.DataFrame) : pandas DataFrame containing the results of
                                        the GRAFIMO analysis
    """

    if not isinstance(sequence_loc, str):
        errmsg = ''.join("\n\nERROR: unable to locate extracted sequences in ", sequence_loc, ". Exiting")
        raise FileNotFoundError(errmsg)

    if not isinstance(motif, Motif):
        raise ValueError("\n\nERROR: the given motif is not an instance of Motif")

    if not isinstance(args_obj, Findmotif):
        raise ValueError("\n\nERROR: unrecognized argument object type")

    # reading arguments
    cores = args_obj.get_cores()
    threshold = args_obj.get_threshold()
    no_qvalue = args_obj.get_no_qvalue()
    qval_t = args_obj.get_qvalueT()
    no_reverse = args_obj.get_no_reverse()
    verbose = args_obj.get_verbose()

    assert threshold > 0
    assert threshold <= 1
    assert cores >= 1

    print_scoring_msg(no_reverse, motif)

    cwd = os.getcwd()
    os.chdir(sequence_loc)  # go to sequence location

    manager = mp.Manager()
    return_dict = manager.dict()  # results
    scanned_nucs_dict = manager.dict()  # nucleotides scanned
    scanned_seqs_dict = manager.dict()  # sequences scanned

    sequences = glob.glob('*.tsv')  # get all tmp files containing sequences
    sequences_split = np.array_split(sequences, cores)  # split the sequence set in #cores chunks

    jobs = []  # jobs list
    proc_finished = 0  # number of jobs done

    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)  # overwrite the default SIGINT handler to exit gracefully
    # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python

    if verbose:
        start_s = time.time()

    try:

        # compute results in parallel
        for i in range(cores):
            p = mp.Process(target=score_seqs, args=(sequences_split[i], motif, no_reverse, return_dict,
                                                         scanned_seqs_dict, scanned_nucs_dict, i))
            jobs.append(p)
            p.start()  # start the process
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
            end_s = time.time()
            msg = ''.join(["\nScored all sequences in ", str(end_s - start_s), "s"])
            print(msg)

        else:
            pass # all was OK, go to the next instruction
        # end if
    # end try

    os.chdir(cwd)  # get back to starting point

    cmd = "rm -rf {0}".format(sequence_loc)  # remove temporary sequence files
    code = subprocess.call(cmd, shell=True)

    if code != 0:
        msg = ' '.join(["\n\nERROR: an error occurred while running", cmd])
        raise SubprocessError(msg)
    # end if

    if verbose:
        start_df = time.time()

    # recover all analysis results and summarize them in a single data-structure
    seqnames = []
    seqs = []
    chroms = []
    starts = []
    stops = []
    strands = []
    scores = []
    pvalues = []
    references = []

    seqs_scanned = 0
    nucs_scanned = 0

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
        references += return_dict[key].get_references()

        # compute the total number of scanned sequences and nucleotides
        seqs_scanned += scanned_seqs_dict[key]  # the keys are the same as return_dict
        nucs_scanned += scanned_nucs_dict[key]  # the keys are the same as return_dict
    # end for

    # compute the q-values
    if no_qvalue:
        qvalues = []  # empty list -> not computed
    else:
        qvalues = compute_qvalues(pvalues)
    # end if

    print("Scanned sequences:", seqs_scanned)
    print("Scanned nucleotides:", nucs_scanned)

    # summarize results in a pandas DF
    finaldf = build_df(motif, seqnames, starts, stops, strands, scores,
                       pvalues, qvalues, seqs, references, threshold,
                       qval_t, no_qvalue)

    if verbose:
        end_df = time.time()
        msg = ''.join(["\nBuilt result summary in ", str(end_df - start_df), "s"])

    return finaldf

# end of compute_results()


def score_seqs(sequences,
               motif,
               no_reverse,
               return_dict,
               scanned_seqs_dict,
               scanned_nucs_dict,
               pid):
    """
        Retrieve the extracted sequences and score them,
        using a given motif
        ----
        Parameters:
            sequences (str) : sequences to score
            motif (Motif) : processed motif used to score sequences
            no_reverse (bool) : if set to True, only sequences belonging
                                to fw strand will be scored, belonging to
                                both fw and rev strand, otherwise
            return_dict (dict) : shared variable used to store partial results
            scanned_seqs_dict (dict) : shared variable to store the number of
                                        scanned sequences
            scanned_nucs_dict (dict) :  shared variable to store the number of
                                        scanned nucleotides
            pid (int) : current process ID
        ----
        Returns:
            None
    """

    try:
        # get Motif attributes to score sequences
        score_matrix = motif.getMotif_scoreMatrix()
        pval_mat = motif.getMotif_pval_mat()
        min_score = motif.getMin_val()
        scale = motif.getScale()
        width = motif.getWidth()
        offset = motif.getOffset()

        # initialize lists where results will be stored
        seqs = []
        scores = []
        pvalues = []
        seqnames = []
        chroms = []
        starts = []
        stops = []
        strands = []
        references = []

        seqs_scanned = 0  # counter of the scanned sequences number

        # read the sequences
        for s in sequences:
            with open(s, mode='r') as raw_sequences:  # read sequence files in read only mode
                for line in raw_sequences:
                    data = line.split('\t')  # are TSV files

                    strand = data[2][-1]

                    if no_reverse:  # score only the fw strand

                        if strand == '+':

                            # read the final values
                            seq = data[1]
                            seqname = data[0]
                            chrom = seqname.split(':')[0]
                            start = data[2].split(':')[1]
                            start = start[:-1]
                            stop = data[3].split(':')[1]
                            stop = stop[:-1]
                            ref = data[4]
                            score, pvalue = compute_score_seq(seq, score_matrix, pval_mat,
                                                              min_score, scale, width, offset)
                            seqs_scanned += 1
                        # end if

                    else:  # score both fw and reverse strands

                        seq = data[1]
                        seqname = data[0]
                        chrom = seqname.split(':')[0]
                        start = data[2].split(':')[1]
                        start = start[:-1]
                        stop = data[3].split(':')[1]
                        stop = stop[:-1]
                        ref = data[4]
                        score, pvalue = compute_score_seq(seq, score_matrix, pval_mat,
                                                          min_score, scale, width, offset)
                        seqs_scanned += 1
                    # end if

                    # add data to final summary
                    seqs.append(seq)
                    scores.append(score)
                    pvalues.append(pvalue)
                    seqnames.append(seqname)
                    chroms.append(chrom)
                    starts.append(start)
                    stops.append(stop)
                    strands.append(strand)
                    references.append(ref)
                # end for
            # end open
        # end for

    except KeyboardInterrupt:
        pass

    else:

        res_tmp = ResultTmp(seqnames, seqs, chroms, starts,
                            stops, strands, scores, pvalues, references)

        return_dict[pid] = res_tmp
        scanned_seqs_dict[pid] = seqs_scanned
        scanned_nucs_dict[pid] = seqs_scanned * width  # width nucleotides per sequence

    # end try

# end of score_seqs()


# using numba are achieved better perfomances
@jit(nopython=True)
def compute_score_seq(seq,
                      score_matrix,
                      pval_mat,
                      min_score,
                      scale,
                      width,
                      offset):
    """
        Score a sequence using a processed motif scoring matrix
        ----
        Parameters:
            seq (str) : sequence to score
            score_matrix (np.ndarray) : scoring matrix
            pval_mat (np.ndarray) : matrix used to compute P-values
                                    using a DP-algorithm (Staden, 1994)
            min_score (int) : lowest score in the scoring matrix
            scale (int) : scale used during motif processing
            width (int) : motif width
            offset (int) : offset used during motif processing
        ----
        Returns:
            score (np.double) : sequence score
            pvalue (np.double) : sequence score P-value
    """

    score = 0

    seq_len = len(seq)
    assert seq_len == width

    # score the current sequence
    for i in range(width):
        nuc = seq[i]

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
    pvalue = (pval_mat[score:].sum()) / tot

    # retrieve the log-likelihood score
    logodds = (score / scale) + (width * offset)
    score = logodds

    assert pvalue > 0
    assert pvalue <= 1

    return score, pvalue

# end of compute_score_seq()


def compute_qvalues(pvalues):
    """
        Compute the q-values for a given list
        of P-values, using the Benjamini-Hochberg method
        ----
        Parameters:
            pvalues (list) : list of P-values
        ----
        Returns:
            qvalues (list) : list of computed q-values
    """

    if not isinstance(pvalues, list):
        errmsg = "\n\nERROR: P-values must be in a list"
        raise ValueException(errmsg)

    print("\nComputing q-values...\n")

    # use Benjamini-Hochberg procedure to correct P-values
    mt_obj = multipletests(pvalues, method="fdr_bh")
    qvalues = list(mt_obj[1])

    return qvalues

# end of compute_qvalues()


def build_df(motif,
             seqnames,
             starts,
             stops,
             strands,
             scores,
             pvalues,
             qvalues,
             sequences,
             references,
             threshold,
             qval_t,
             no_qvalue):
    """
        Build a pandas DataFrame to summarize the results
        of GRAFIMO analysis
        ----
        Parameters:
            motif (Motif) : motif
            seqnames (list) : list of sequence names
            starts (list) : list of sequence starting positions
            stops (list) : list of sequence ending positions
            strands (list) : list of sequence strands
            scores (list) : list of sequence scores
            pvalues (list) : list of sequence score P-values
            qvalues (list) : list of sequence q-values
            sequences (list) : list of sequences
            references (list) : list of sequence flag values. If 'ref',
                                then the sequence belong to the reference genome,
                                if 'non.ref', then the sequence contains variants
            threshold (float) : threshold to apply on the P-value (default behavior)
                                or on the q-values
            qval_t (bool) : if set to True, the threshold will be applied on the
                            q-values, on the P-values otherwise
        ----
        Returns:
             df (pd.DataFrame)
    """

    if not isinstance(motif, Motif):
        errmsg = "\n\nERROR: unknown data-type for motif"
        raise ValueException(errmsg)

    if not isinstance(seqnames, list):
        errmsg = "\n\nERROR: unknown data-type, cannot proceed"
        raise ValueException(errmsg)

    if not isinstance(starts, list):
        errmsg = "\n\nERROR: unknown data-type, cannot proceed"
        raise ValueException(errmsg)

    if not isinstance(stops, list):
        errmsg = "\n\nERROR: unknown data-type, cannot proceed"
        raise ValueException(errmsg)

    if not isinstance(strands, list):
        errmsg = "\n\nERROR: unknown data-type, cannot proceed"
        raise ValueException(errmsg)

    if not isinstance(pvalues, list):
        errmsg = "\n\nERROR: unknown data-type, cannot proceed"
        raise ValueException(errmsg)

    if not isinstance(qvalues, list):
        errmsg = "\n\nERROR: unknown data-type, cannot proceed"
        raise ValueException(errmsg)

    if not isinstance(sequences, list):
        errmsg = "\n\nERROR: unknown data-type, cannot proceed"
        raise ValueException(errmsg)

    if not isinstance(references, list):
        errmsg = "\n\nERROR: unknown data-type, cannot proceed"
        raise ValueException(errmsg)

    if not isinstance(qval_t, bool):
        errmsg = "\n\nERROR: unknown data-type, cannot proceed"
        raise ValueException(errmsg)

    if not isinstance(no_qvalue, bool):
        errmsg = "\n\nERROR: unknown data-type, cannot proceed"
        raise ValueException(errmsg)

    # all lists must have the same length
    lst_len = len(seqnames)

    assert len(starts) == lst_len
    assert len(stops) == lst_len
    assert len(strands) == lst_len
    assert len(scores) == lst_len
    assert len(pvalues) == lst_len
    assert len(sequences) == lst_len
    assert len(references) == lst_len

    # check if we want also the q-values
    if not no_qvalue:  # we want the q-values
        assert len(qvalues) == lst_len

    if qval_t:  # apply the threshold on the q-values rather than on P-values
        assert (not no_qvalue)
        assert len(qvalues) > 0   # we must have computed them

    seqnames_thresh = []
    starts_thresh = []
    ends_thresh = []
    strands_thresh = []
    scores_thresh = []
    pvalues_thresh = []
    sequences_thresh = []
    references_thresh = []

    if not no_qvalue:
        qvalues_thresh = []

    for i in range(lst_len):

        if not qval_t:  # apply threshold on P-values
            pvalue = pvalues[i]

            if pvalue < threshold:

                # only the sequences with a P-value under the threshold survive
                seqnames_thresh.append(seqnames[i])
                starts_thresh.append(starts[i])
                ends_thresh.append(stops[i])
                strands_thresh.append(strands[i])
                scores_thresh.append(scores[i])
                pvalues_thresh.append(pvalues[i])
                sequences_thresh.append(sequences[i])
                references_thresh.append(references[i])

                if not no_qvalue:
                    qvalues_thresh.append(qvalues[i])
            # end if

        else:  # apply threshold on q-values
            qvalue = qvalues[i]

            if qvalue < threshold:

                # only the sequences with a q-value under the threshold survive
                seqnames_thresh.append(seqnames[i])
                starts_thresh.append(starts[i])
                ends_thresh.append(stops[i])
                strands_thresh.append(strands[i])
                scores_thresh.append(scores[i])
                pvalues_thresh.append(pvalues[i])
                sequences_thresh.append(sequences[i])
                references_thresh.append(references[i])

                # the last control statement, in the if, in this case is not
                # necessary (we must have the q-values)
                # otherwise we should not be here
                qvalues_thresh.append(qvalues[i])
            # end if
        # end if
    # end for

    df_len = len(seqnames_thresh)

    # TF's name and ID list
    motif_ids = [motif.getMotifID()] * df_len
    motif_names = [motif.getMotifName()] * df_len

    """
        build the final data frame
        
        structure:
        
           |motif_id|motif_alt_id|sequence_name|start|stop|strand|score|p-value|q-value|matched_sequence|reference|   
    """

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
    df['reference'] = references_thresh

    # sort entries by p-value
    df = df.sort_values(['p-value'], ascending=True)

    # reindex the data frame in order to have indexes in range [1, (df_len + 1)]
    df.index = list(range(1, (df_len + 1)))

    return df

# end of build_df()


def print_scoring_msg(no_reverse, motif):
    if not isinstance(motif, Motif):
        raise ValueException('\n\nERROR: The given motif is not an instance of Motif')

    motif_id = motif.getMotifID()
    fw_id = ''.join(['+', motif_id])

    # we take into account also the reverse complement
    if not no_reverse:
        rev_id = ''.join(['-', motif_id])

    print()  # newline
    print('Scoring hits for motif', fw_id)

    # if we score also the reverse complement
    if not no_reverse:
        print('Scoring hits for motif', rev_id)

    print()  # newline

# end of print_scoring_msg()

