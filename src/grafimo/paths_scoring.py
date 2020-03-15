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

import pandas as pd
from . import motif as mtf
from grafimo.GRAFIMOException import WrongPathException, ValueException, SubprocessException
from grafimo.utils import die, printProgressBar
import os
import multiprocessing as mp
import numpy as np
import glob
import subprocess
from numba import jit
from statsmodels.stats.multitest import multipletests

def scoreGraphsPaths(subgraphs, motif, threshold, cores, no_reverse, qvalue, qvalueT):
    """
        Score the extracted sequences and save the results in a pandas 
        DataFrame. 
        ----
        Parameters:
            subgraphs (str) : path to the subgraphs sequences
            motif (Motif) :
            threshold (float) : threshold to apply on the p-value (default) or on
                                the q-value
            cores (int) : cores to use to run the function in parallel
            no_reverse (bool) : flag to consider only the sequences from the
                                forward strand
            qvalue (bool) : flag to compute also the q-value for each hit
                            found
            qvalueT (bool) : flag value, if set to True, the threshold will be
                                applied on the q-values
        ----
        Returns:
            finaldf (pd.DataFrame) : dataframes with the results
            
    """
    
    if not isinstance(subgraphs, str):
        raise WrongPathException("The given path to the subgraphs is not valid")
        die(1)
    
    if not isinstance(motif, mtf.Motif):
        raise ValueException("The given motif is not an instance of Motif")
        die(1)

    if not isinstance(no_reverse, bool):
        raise ValueException("Do not know what to do with this 'no_reverse' value")
        die(1)

    if not isinstance(qvalue, bool):
        raise ValueException("Do not know what to do with this 'qvalue' value")
        die(1)
        
    assert threshold > 0
    assert threshold <= 1
    assert cores >= 0

    # the args are correct
        
    scoringPathsMsg(no_reverse, motif)
        
    if cores == 0:
        NCORES = mp.cpu_count()
    else:
        NCORES = cores

    assert NCORES > 0
        
    cwd = os.getcwd()
    os.chdir(subgraphs)
        
    manager = mp.Manager()
    returnDict = manager.dict() # results
    scannedNucsDict = manager.dict() # nucleotides scanned
    scannedSeqsDict = manager.dict() # sequences scanned
        
    sgs = glob.glob('*.tsv')
    sgs_splt = np.array_split(sgs, NCORES)
    jobs = []
    proc_finished = 0

    # run processes in parallel
    for i in range(NCORES):
        p = mp.Process(target=score_subgraphs, args=(sgs_splt[i], motif, no_reverse, i, returnDict,
                                                        scannedSeqsDict, scannedNucsDict))
        jobs.append(p)
        p.start() # start the process
            
    printProgressBar(proc_finished, NCORES, prefix='Progress:',
                            suffix='Complete', length=50)
    for job in jobs:
        job.join() # deadlock
        proc_finished += 1
        printProgressBar(proc_finished, NCORES, prefix='Progress:',
                            suffix='Complete', length=50)
            
    os.chdir(cwd)
        
    cmd = "rm -rf {0}".format(subgraphs)
    code = subprocess.call(cmd, shell=True)
        
    if code != 0:
        msg = ' '.join(["The command", cmd, "returned non-zero exit status"])
        raise SubprocessException(msg)
        die(1)

    # initialize the dataframe columns
    seqs = []
    scores = []
    seqnames = []
    chroms = []
    starts = []
    ends = []
    pvalues = []
    qvalues = []
    strands = []
    references = []

    for key in returnDict.keys():

        seqs += returnDict[key]['sequence'].to_list()
        scores += returnDict[key]['score'].to_list()
        seqnames += returnDict[key]['sequence_name'].to_list()
        chroms += returnDict[key]['chromosome'].to_list()
        starts += returnDict[key]['start'].to_list()
        ends += returnDict[key]['end'].to_list()
        pvalues += returnDict[key]['pvalue'].to_list()
        strands += returnDict[key]['strand'].to_list()
        references += returnDict[key]['reference'].to_list()

    # get the total sequences scanned
    seqsScanned = 0
    for key in scannedSeqsDict.keys():
        seqsScanned += scannedSeqsDict[key]

    # get the total nucleotides scanned
    nucsScanned = 0
    for key in scannedNucsDict.keys():
        nucsScanned += scannedNucsDict[key]

    # compute the q-values
    if qvalue:
        qvalues = compute_qvalues(pvalues)

    print("\nScanned sequences:", seqsScanned)
    print("Scanned nucleotides:", nucsScanned)
            
    finaldf = buildDF(motif, seqnames, starts, ends, strands, scores,
                        pvalues, qvalues, seqs, references, threshold, qvalueT)
            
    return finaldf


def get_subgraphs_dict(sgs):
    """
        Keep track of the strand to which each sequence belong
        ----
        Parameters:
            sgs (str) : path to the subgraphs
        ----
        Returns:
            sg_seqs_dict (dict) : dictionary of all the sequences
            sg_directions_dict (dict) : direction of each sequence
            sg_starts_dict (dict) : starting positions of the sequences
            sg_ends_dict (dict) : ending positions of the sequences
    """
    
    sg_seqs_dict = {} # sequences
    sg_directions_dict = {} # direction (forward/reverse)
    sg_starts_dict = {} # starting positions
    sg_ends_dict = {} # ending positions
   
    for sg in sgs:
        sg_name = sg.split('.')[0]
        sg_data = pd.read_csv(sg, header=None, sep='\t')
        sg_seqs_dict.update({sg_name : sg_data.loc[:,1].to_list()})
        
        seq_directions = get_paths_directions(sg_data)
        sg_directions_dict.update({sg_name : seq_directions})

        starts, ends = get_paths_positions(sg_data)
        sg_starts_dict.update({sg_name : starts})
        sg_ends_dict.update({sg_name : ends})
        
    return sg_seqs_dict, sg_directions_dict, \
                sg_starts_dict, sg_ends_dict


def get_paths_directions(paths):
    
    paths_num=len(paths.index)
    directions=allocate_array(None, paths_num)
    seq_dir_list=paths.loc[:,2].to_list() # retrieve the direction field
    
    
    for i in range(paths_num):
        seq_dir = seq_dir_list[i]
        
        if '-' in seq_dir:
            directions[i] = 1 # reverse strand
            
        else:
            directions[i] = 0 # forward strand
            
    return directions


def get_paths_positions(paths):

    paths_num = len(paths.index)
    starts = allocate_array(None, paths_num)
    ends = allocate_array(None, paths_num)
    seq_start_list = paths.loc[:, 2].to_list() # retrieve starts
    seq_end_list = paths.loc[:, 3].to_list() # retrieve ends

    for i in range(paths_num):
        # get start
        seq_start = seq_start_list[i].split(':')[1] 
        seq_start = seq_start[:-1]
        starts[i] = seq_start

        #get end
        seq_end = seq_end_list[i].split(':')[1]
        seq_end = seq_end[:-1]
        ends[i] = seq_end

    return starts, ends


def allocate_array(value=None, size='0'):
    
    arr = [value] * size
    
    return arr
    
        
def score_subgraphs(sgs, motif, no_reverse, psid, returnDict,
                        scannedSeqsDict, scannedNucsDict):
    """
        Score all the sequences contained in the extracted
        subgraphs.
        ----
        Parameters:
            sgs (str) : path to the subgraphs
            motif (Motif) : Motif object
            no_reverse (bool) : flag to consider only the forward strand or
                                both the forard and reverse strands
            psid (int) : process ID of the scoring process, defined during the
                            parallelization step
            returnDict (dict) : dictionary to reconstruct the output after the
                                parallelization
            scannedSeqsDict (dict) : dictionary to keep track of the sequences scanned
            scannedNucsDict (dict) : dictionary to keep track of the nucleotides scanned
        ----
        Returns:
            None
    """

    # get attributes to score the sequences
    scoreMatrix = motif.getMotif_scoreMatrix()
    pval_mat = motif.getMotif_pval_mat()
    minScore = motif.getMin_val()
    scale = motif.getScale()
    width = motif.getWidth()
    offset = motif.getOffset()

    # lists for the final results 
    seqs = []
    scores_final = []
    pvalues_final = []
    seqnames = []
    chroms = []
    starts = []
    ends = []
    strands = []
    references = []

    seqsScanned = 0 # counter for the sequences scanned
    
    for sg in sgs:
        
        # get the sequences and their starting positions of the  
        # current subgraph
        sg_data = pd.read_csv(sg, header = None, sep = '\t')
        seqs_sg = np.asarray(sg_data[1].to_list())
        starts_sg = np.asarray(sg_data[2].to_list())

        seqs_num = len(seqs_sg)

        assert seqs_num == len(starts_sg)
        
        # score the sequences
        scores, pvalues = score_sequences(seqs_sg, starts_sg, no_reverse, scoreMatrix, 
                                            pval_mat, minScore, scale, width, offset, seqs_num)
        scores_final += scores
        pvalues_final += pvalues
        
        # get the remaining informations
        for i in range(seqs_num):

            strand = starts_sg[i][-1]

            if no_reverse:
                # we are interested only in the sequences coming form the forward strand
                if strand == '+': # is forward
                    seq = sg_data.loc[i, 1]
                    seqname = sg_data.loc[i, 0]
                    chrom = seqname.split(':')[0]
                    start = sg_data.loc[i, 2].split(':')[1]
                    start = start[:-1]
                    end = sg_data.loc[i, 3].split(':')[1]
                    end = end[:-1]
                    reference = sg_data.loc[i, 4]

                    seqsScanned += 1

                    seqs.append(seq)
                    seqnames.append(seqname)
                    strands.append(strand)
                    chroms.append(chrom)
                    starts.append(start)
                    ends.append(end)
                    references.append(reference)
                        
            else:
                # consider both strands
                seq = sg_data.loc[i, 1]
                seqname = sg_data.loc[i, 0]
                chrom = seqname.split(':')[0]
                start = sg_data.loc[i, 2].split(':')[1]
                start = start[:-1]
                end = sg_data.loc[i, 3].split(':')[1]
                end = end[:-1]
                reference = sg_data.loc[i, 4]

                seqsScanned += 1

                seqs.append(seq)
                seqnames.append(seqname)
                strands.append(strand)
                chroms.append(chrom)
                starts.append(start)
                ends.append(end)
                references.append(reference)


    # build the dataframe for the results got from the current process
    summary = pd.DataFrame()
    summary['sequence_name'] = seqnames
    summary['chromosome'] = chroms
    summary['start'] = starts
    summary['end'] = ends
    summary['sequence'] = seqs
    summary['score'] = scores_final
    summary['pvalue'] = pvalues_final
    summary['strand'] = strands
    summary['reference'] = references
    
    # retrieve the final results summary
    returnDict[psid] = summary
    scannedSeqsDict[psid] = seqsScanned
    scannedNucsDict[psid] = seqsScanned * motif.getWidth() # in every sequence we have scanned motif_width nucleotides


# here is used numba to achieve performances near to the ones 
# obtained with a C program
@jit(nopython=True)
def score_sequences(seqs, starts, no_reverse, scoreMatrix, 
                        pval_mat, minScore, scale, width, offset, seqs_num):
    """
        Scores each sequence contained in the list seqs, using the values 
        defined in the score matrix, which were obtained from the 
        input motif file.
        ----
        Parameters:
            seqs (list) : list of DNA sequences to score
            starts (list) : list of the starting positions 
                            of the sequences contained in seqs
            no_reverse (bool) : flag value to consider only
                                sequences from the forward strand
            scoreMatrix (np.ndarray) : scoring matrix 
            pval_mat (np.array) : p-value matrix
            minScore (np.double) : minimum score of the scoring
                                    matrix
            scale (int) : scale used during the scaling of the
                            scoring matrix
            width (int) : width of the motif
            offset (np.double) : offset to retrieve the log 
                                    likelihood score from the
                                    scaled score
            seqs_num (int) : number of the sequences to score
        ----
        Returns:
            scores (list) : list of the scores of the sequences
            pvalues (list) : list of the p-values of the 
                                sequences  
    """

    scores = []
    pvalues = []

    # loop over every sequnce    
    for i in range(seqs_num):

        strand = starts[i][len(starts[i]) - 1] # get the strand

        if no_reverse:
            if strand == '+': # is forward
                seq = seqs[i]
                score = 0

                seq_len = len(seq)

                # score the current sequence
                for i in range(seq_len):
                    nuc = seq[i]

                    if nuc == 'N':
                        score = minScore
                        break # we don't go further

                    if nuc == 'A' or nuc == 'a':
                        score += scoreMatrix[0, i]
                    if nuc == 'C' or nuc == 'c':
                        score += scoreMatrix[1, i]
                    if nuc == 'G' or nuc == 'g':
                        score += scoreMatrix[2, i]
                    if nuc == 'T' or nuc == 't':
                        score += scoreMatrix[3, i]

                # get the p-value for the obtained score
                tot = pval_mat.sum()
                pvalue = (pval_mat[score:].sum())/tot

                # retrieve the loglikelihood score
                logodds = (score / scale) + (width * offset)
                score = logodds
                
                scores.append(score)
                pvalues.append(pvalue)
        
        else:

            seq = seqs[i]
            score = 0

            seq_len = len(seq)

            # score the current sequence
            for i in range(seq_len):
                nuc = seq[i]

                if nuc == 'N':
                    score = minScore
                    break # we don't go further

                if nuc == 'A' or nuc == 'a':
                    score += scoreMatrix[0, i]
                if nuc == 'C' or nuc == 'c':
                    score += scoreMatrix[1, i]
                if nuc == 'G' or nuc == 'g':
                    score += scoreMatrix[2, i]
                if nuc == 'T' or nuc == 't':
                    score += scoreMatrix[3, i]

            # get the p-value for the obtained score
            tot = pval_mat.sum()
            pvalue = (pval_mat[score:].sum())/tot

            # retrieve the loglikelihood score
            logodds = (score / scale) + (width * offset)
            score = logodds
                
            scores.append(score)
            pvalues.append(pvalue)
        
        # end if
    # end for

    return scores, pvalues


def compute_qvalues(pvalues):
    """
        Compute the corresponding q-values for a given list
        of p-values, using the Benjamini-Hockberg method
        ----
        Parameters:
            pvalues (list) : list of p-values
        ----
        Returns:
            qvalues (list) : list of computed q-values
    """

    if not isinstance(pvalues, list):
        raise ValueException("To compute the corresponding p-values, the q-values must be stored in a list")
        die(1)

    print("\nComputing q-values...\n")

    mt_obj = multipletests(pvalues, method="fdr_bh")
    qvalues = list(mt_obj[1])

    return qvalues

        
def buildDF(motif, seqnames, starts, ends, strands, 
                scores, pvalues, qvalues, sequences,
                references, threshold, qvalueT):
    """
        Build the dataframe summary of the retrieved data
        ----
        Params:
            motif (Motif) : Motif object
            seqnames (list) : list of the sequence names
            starts (list) : list of the starting positions of the hits
            ends (list) : list of the ending positions of the hits
            strands (list) : list of the strands of the hits
            scores (list) : list of the scores of the hits
            pvalues (list) : list of the p-values for the hits
            qvalues (list) : list of q-values for the hits
            sequences (list) : list of the hits
            references (list) : list that keep track of the fact if a sequence
                                belongs to the reference genome or is obtained
                                from the variants dfined in the VCF used to
                                build the genome graph
            threshold (float) : threshold to apply on the p-value (default) or on
                                the q-values
            qvalueT (bool) : if set to True, the threshold will be applied on the
                                q-values
        ----
        Returns:
             df (pd.DataFrame)
    """
    
    if not isinstance(motif, mtf.Motif):
        raise ValueException("The motif is not an instance of Motif")
        die(1)
        
    if not isinstance(seqnames, list):
        raise ValueException("The sequence names must be in a list")
        die(1)
        
    if not isinstance(starts, list):
        raise ValueException("The starts must be in a list")
        die(1)
        
    if not isinstance(ends, list):
        raise ValueException("The ends must be in a list")
        die(1)
        
    if not isinstance(strands, list):
        raise ValueException("The strands must be in a list")
        die(1)
        
    if not isinstance(pvalues, list):
        raise ValueException("The p-values must be in a list")
        die(1)

    if not isinstance(qvalues, list):
        raise ValueException("The q-values must be in a list")
        die(1)
        
    if not isinstance(sequences, list):
        raise ValueException("The sequences must be in a list")
        die(1)

    if not isinstance(references, list):
        raise ValueException("The reference list must be of list type")

    # all lists must have the same length
    LSTLEN = len(seqnames)

    assert len(starts) == LSTLEN
    assert len(ends) == LSTLEN
    assert len(strands) == LSTLEN
    assert len(scores) == LSTLEN
    assert len(pvalues) == LSTLEN
    assert len(sequences) == LSTLEN
    assert len(references) == LSTLEN

    # check if we want also the q-values
    QVAL = False
    if len(qvalues) > 0: # we want the q-values
        QVAL = True
        assert len(qvalues) == LSTLEN

    if qvalueT:
        assert len(qvalues) > 0 # if we apply the threshold on the q-values
                                # we must have computed them

    seqnames_thresh = []
    starts_thresh = []
    ends_thresh = []
    strands_thresh = []
    scores_thresh = []
    pvalues_thresh = []
    sequences_thresh = []
    references_thresh = []

    if QVAL:
        qvalues_thresh = []

    for idx in range(LSTLEN):

        if not qvalueT:
            pvalue = pvalues[idx]
            if pvalue < threshold:
                # only the sequences with a p-value under the threshold survive
                seqnames_thresh.append(seqnames[idx])
                starts_thresh.append(starts[idx])
                ends_thresh.append(ends[idx])
                strands_thresh.append(strands[idx])
                scores_thresh.append(scores[idx])
                pvalues_thresh.append(pvalues[idx])
                sequences_thresh.append(sequences[idx])
                references_thresh.append(references[idx])

                if QVAL:
                    qvalues_thresh.append(qvalues[idx])

        else:
            qvalue = qvalues[idx]
            if qvalue < threshold:
                # only the sequences with a q-value under the threshold survive
                seqnames_thresh.append(seqnames[idx])
                starts_thresh.append(starts[idx])
                ends_thresh.append(ends[idx])
                strands_thresh.append(strands[idx])
                scores_thresh.append(scores[idx])
                pvalues_thresh.append(pvalues[idx])
                sequences_thresh.append(sequences[idx])
                references_thresh.append(references[idx])
                qvalues_thresh.append(qvalues[idx]) # the last control statement in this case is not
                                                    # necessary (we must have the q-values)
                                                    # otherwise we should not be here

    DFLEN = len(seqnames_thresh)

    # TF's name and ID list
    motifIDs = [motif.getMotifID()] * DFLEN
    motifNames = [motif.getMotifName()] * DFLEN

    # build the final dataframe
    df = pd.DataFrame()
    df['motif_id'] = motifIDs
    df['motif_alt_id'] = motifNames
    df['sequence_name'] = seqnames_thresh
    df['start'] = starts_thresh
    df['stop'] = ends_thresh
    df['strand'] = strands_thresh
    df['score'] = scores_thresh
    df['p-value'] = pvalues_thresh

    # add the q-values to the final dataframe if they have been computed
    if QVAL:
        df['q-value'] = qvalues_thresh

    # finish to build the data frame
    df['matched_sequence'] = sequences_thresh
    df['reference'] = references_thresh
        
    # values are sorted by p-value
    df = df.sort_values(['p-value'], ascending=True)

    # reindex the dataframe in order to have indexes from 1 to DFLEN + 1
    df.index = list(range(1, DFLEN + 1))
        
    return df
    

def scoringPathsMsg(no_reverse, motif):
    
    if not isinstance(motif, mtf.Motif):
        raise ValueException('The given motif is not an instance of Motif. Cannot get its ID')
        die(1)
        
    motifID=motif.getMotifID()
    posID=''.join(['+', motifID])

    # we take into account also the reverse complement
    if not no_reverse:
        revID=''.join(['-', motifID])

    print() # newline
    print('Scoring hits for motif', posID, sep=' ')
    
    # if we score also the reverse complement
    if not no_reverse: 
        print('Scoring hits for motif', revID, sep=' ')
    
    print() # newline
    
