"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Script that computes the scores, for each extracted sequence from the peak 
defined in the bedfile and obtained from the genome graph, and saves the 
results in a pandas DataFrame

"""

import pandas as pd
from . import motif as mtf
from grafimo.GRAFIMOException import WrongPathException, ValueException, SubprocessException
from grafimo.utils import die, printProgressBar
from GRAFIMOscoring import score_seq
import os
import multiprocessing as mp
import numpy as np
import glob
import subprocess
from statsmodels.stats.multitest import multipletests

def scoreGraphsPaths(subgraphs, motif, pvalueT, cores, no_reverse, qvalue):
    """
        Score the extracted sequences and save the results in a pandas 
        DataFrame. A parallel computing approach is used.
        ----
        Parameters:
            subgraphs (str) : path to the subgraphs sequences
            motif (Motif) :
            pvalueT (float) : p-value threshold
            cores (int) : cores to use to run the function in parallel
            no_reverse (bool) : flag to consider only the sequences from the
                                forward strand
            qvalue (bool) : flag to compute also the q-value for each hit
                            found
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
        
    assert pvalueT > 0
    assert cores >= 0

    # the args are correct
        
    scoringPathsMsg(no_reverse, motif)
        
    if cores==0:
        NCORES=mp.cpu_count()
    else:
        NCORES=cores

    assert NCORES > 0
        
    cwd=os.getcwd()
    os.chdir(subgraphs)
        
    manager = mp.Manager()
    returnDict = manager.dict() # results
    scannedNucsDict = manager.dict() # nucleotides scanned
    scannedSeqsDict = manager.dict() # sequences scanned
        
    sgs=glob.glob('*.tsv')
    sgs_splt=np.array_split(sgs, NCORES)
    jobs=[]
    proc_finished=0

    # run the processes in parallel
    for i in range(NCORES):
        p=mp.Process(target=score_subgraphs, args=(sgs_splt[i], motif, no_reverse, i, returnDict,
                                                    scannedSeqsDict, scannedNucsDict))
        jobs.append(p)
        p.start() # start the process
            
    for job in jobs:
        proc_finished += 1
        printProgressBar(proc_finished, NCORES, prefix='Progress:',
                            suffix='Complete', length=100)
        job.join() # deadlock
            
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
        qvalues = []
        qvalues = compute_qvalues(pvalues)

    print("\nScanned sequences:", seqsScanned)
    print("Scanned nucleotides:", nucsScanned)
            
    finaldf = buildDF(motif, seqnames, starts, ends, strands, scores,
                        pvalues, qvalues, seqs, references, pvalueT)
            
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
        Score the sequences extracted from the subgraphs
        ----
        Parameters:
            sgs (str) : path to the subgraphs
            motif (Motif) : motif object
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

    scoreMatrix = motif.getMotif_scoreMatrix()
    pval_mat = motif.getMotif_pval_mat()
    minScore = motif.getMin_val()
    scale = motif.getScale()
    width = motif.getWidth()
    offset = motif.getOffset()

    seqs = []
    scores = []
    pvalues = []
    seqnames = []
    chroms = []
    starts = []
    ends = []
    strands = []
    references = []

    seqsScanned = 0 # counter for the sequences scanned
    
    for sg in sgs:
        
        sg_data = pd.read_csv(sg, header = None, sep = '\t')
        paths_num = len(sg_data.index)
        
        for i in range(paths_num):
            strand = sg_data.loc[i, 2][-1]
            if no_reverse:
                if strand == '+': # is forward
                    seq = sg_data.loc[i, 1]
                    seqname = sg_data.loc[i, 0]
                    chrom = seqname.split(':')[0]
                    start = sg_data.loc[i, 2].split(':')[1]
                    start = start[:-1]
                    end = sg_data.loc[i, 3].split(':')[1]
                    end = end[:-1]
                    reference = sg_data.loc[i, 4]
                    score, pvalue = score_seq(seq, scoreMatrix, pval_mat, minScore, scale,
                                                width, offset)
                    # forward strand

                    seqs.append(seq)
                    scores.append(score)
                    pvalues.append(pvalue)
                    seqnames.append(seqname)
                    strands.append(strand)
                    chroms.append(chrom)
                    starts.append(start)
                    ends.append(end)
                    references.append(reference)

                    seqsScanned += 1
                        
            else:

                seq = sg_data.loc[i, 1]
                seqname = sg_data.loc[i, 0]
                chrom = seqname.split(':')[0]
                start = sg_data.at[i, 2].split(':')[1]
                start = start[:-1]
                end = sg_data.loc[i, 3].split(':')[1]
                end = end[:-1]
                reference = sg_data.loc[i, 4]
                score, pvalue = score_seq(seq, scoreMatrix, pval_mat, minScore, scale,
                                            width, offset)
                seqs.append(seq)
                scores.append(score)
                pvalues.append(pvalue)
                seqnames.append(seqname)
                strands.append(strand)
                chroms.append(chrom)
                starts.append(start)
                ends.append(end)
                references.append(reference)

                seqsScanned += 1

                                   
    summary = pd.DataFrame()
    summary['sequence_name'] = seqnames
    summary['chromosome'] = chroms
    summary['start'] = starts
    summary['end'] = ends
    summary['sequence'] = seqs
    summary['score'] = scores
    summary['pvalue'] = pvalues
    summary['strand'] = strands
    summary['reference'] = references
    
    returnDict[psid] = summary
    scannedSeqsDict[psid] = seqsScanned
    scannedNucsDict[psid] = seqsScanned * motif.getWidth() # in every sequence we have scanned motif_width nucleotides


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
                scores, pvalues, qvalues, sequences, references, pvalueT):
    """
        Build the dataframe summarizing the obtained results
        ----
        Params:
            motif (Motif) : motif object
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
            pvalueT (float) : p-value threshold
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
    LSTLEN=len(seqnames)

    assert len(starts) == LSTLEN
    assert len(ends) == LSTLEN
    assert len(strands) == LSTLEN
    assert len(scores) == LSTLEN
    assert len(pvalues) == LSTLEN
    assert len(sequences) == LSTLEN
    assert len(references) == LSTLEN

    # check if we have also the q-values
    QVAL = False
    if len(qvalues) > 0: # there are q-values
        QVAL = True
        assert len(qvalues) == LSTLEN

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
        pvalue = pvalues[idx]

        if pvalue < pvalueT:
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

    DFLEN = len(seqnames_thresh)

    # TF's name and ID list
    motifIDs = [motif.getMotifID()]*DFLEN
    motifNames = [motif.getMotifName()]*DFLEN

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

    df.index = list(range(1, DFLEN+1))
        
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
    
        