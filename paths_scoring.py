"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Script that computes the scores, for each extracted sequence from the peak 
defined in the bedfile and obtained from the genome graph, and saves the 
results in a pandas DataFrame

"""

#### add scanned sequences
#### add scanned positions

import pandas as pd
from . import motif as mtf
from grafimo.GRAFIMOException import WrongPathException, ValueException, SubprocessException
from grafimo.utils import die
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
        N_CORES=mp.cpu_count()
    else:
        N_CORES=cores

    assert N_CORES > 0
        
    cwd=os.getcwd()
    os.chdir(subgraphs)
        
    manager = mp.Manager()
    returnDict = manager.dict() # results
    seenNucsDict = manager.dict() # nucleotides seen
    seenSeqsDict = manager.dict() # sequences seen
        
    sgs=glob.glob('*.tsv')
    sgs_splt=np.array_split(sgs, N_CORES)
    jobs=[]
    proc_finished=0

    # run the processes in parallel
    for i in range(N_CORES):
        p=mp.Process(target=score_subgraphs, args=(sgs_splt[i], motif, no_reverse, i, returnDict,
                                                    seenSeqsDict, seenNucsDict))
        jobs.append(p)
        p.start() # start the process
            
    for job in jobs:
        proc_finished += 1
        printProgressBar(proc_finished, N_CORES, prefix='Progress:',
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

    for key in returnDict.keys():
        seqs += returnDict[key]['sequence'].to_list()
        scores += returnDict[key]['score'].to_list()
        seqnames += returnDict[key]['sequence_name'].to_list()
        chroms += returnDict[key]['chromosome'].to_list()
        starts += returnDict[key]['start'].to_list()
        ends += returnDict[key]['end'].to_list()
        pvalues += returnDict[key]['pvalue'].to_list()
        strands += returnDict[key]['strand'].to_list()

    # get the total sequences seen
    seqsSeen = 0
    for key in seenSeqsDict.keys():
        seqsSeen += seenSeqsDict[key]

    # get the total nucleotides seen
    nucsSeen = 0
    for key in seenNucsDict.keys():
        nucsSeen += seenNucsDict[key]

    # compute the q-values
    qvalues = []
    if qvalue:
        qvalues = compute_qvalues(pvalues)

    print("\nSeen sequences:", seqsSeen)
    print("Seen nucleotides:", nucsSeen)
            
    finaldf = buildDF(motif, seqnames, starts, ends, strands, scores,
                        pvalues, qvalues, seqs, pvalueT)
            
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
    """
    
    sg_seqs_dict={} # sequences
    sg_directions_dict={} # direction (forward/reverse)
   
    for sg in sgs:
        sg_name=sg.split('.')[0]
        sg_data=pd.read_csv(sg, header=None, sep='\t')
        sg_seqs_dict.update({sg_name:sg_data.loc[:,0].to_list()})
        seq_directions=get_paths_directions(sg_data)
            
        sg_directions_dict.update({sg_name:seq_directions})
        
    return sg_seqs_dict, sg_directions_dict


def get_paths_directions(paths):
    
    paths_num=len(paths.index)
    directions=allocate_array(None, paths_num)
    seq_dir_list=paths.loc[:,1].to_list() # retrieve the direction field
    
    
    for i in range(paths_num):
        seq_dir = seq_dir_list[i]
        
        if '-' in seq_dir:
            directions[i] = 1 # reverse strand
            
        else:
            directions[i] = 0 # forward strand
            
    return directions


def allocate_array(value=None, size='0'):
    
    arr = [value] * size
    
    return arr
    
        
def score_subgraphs(sgs, motif, no_reverse, psid, returnDict,
                        seenSeqsDict, seenNucsDict):
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
            seenSeqsDict (dict) : dictionary to keep track of the sequences seen
            seenNucsDict (dict) : dictionary to keep track of the nucleotides seen
        ----
        Returns:
            None
    """
    
    sg_paths_dict, sg_directions_dict = get_subgraphs_dict(sgs)

    # there is no match beetwen the two dictionaries
    if set(sg_paths_dict.keys()) != set(sg_directions_dict.keys()):
        raise ValueException("No match between paths and strands")
        die(1)

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

    seqsSeen = 0 # counter for the sequences seen
    
    for key in sg_paths_dict.keys():
        
        paths = sg_paths_dict[key]
        dirs = sg_directions_dict[key]
        
        for i in range(len(paths)):
            if no_reverse:
                if int(dirs[i]) == 0: # is forward
                    seq = paths[i]
                    seqname = str(key)
                    chrom, pos = seqname.split('_')
                    start, end = pos.split('-')
                    seqname = chrom+':'+start+'-'+end
                    score, pvalue = score_seq(seq, scoreMatrix, pval_mat, minScore, scale,
                                                width, offset)
                    strand= '+' # forward strand

                    seqs.append(seq)
                    scores.append(score)
                    pvalues.append(pvalue)
                    seqnames.append(seqname)
                    strands.append(strand)
                    chroms.append(chrom)
                    starts.append(start)
                    ends.append(end)

                    seqsSeen += 1
                        
            else:
                if int(dirs[i]) == 0: #is forward
                    seq = paths[i]
                    seqname = str(key)
                    chrom, pos = seqname.split('_')
                    start, end = pos.split('-')
                    seqname = chrom+':'+start+'-'+end
                    score, pvalue = score_seq(seq, scoreMatrix, pval_mat, minScore, scale,
                                                width, offset)
                    strand = '+'

                elif int(dirs[i]) == 1: # is reverse
                    seq = paths[i]
                    seqname = str(key)
                    chrom, pos = seqname.split('_')
                    start, end = pos.split('-')
                    seqname = chrom+':' + start + '-'+end
                    score, pvalue = score_seq(seq, scoreMatrix, pval_mat, minScore, scale,
                                                width, offset)
                    strand = '-'

                seqs.append(seq)
                scores.append(score)
                pvalues.append(pvalue)
                seqnames.append(seqname)
                strands.append(strand)
                chroms.append(chrom)
                starts.append(start)
                ends.append(end)

                seqsSeen += 1
                    
    summary = pd.DataFrame()
    summary['sequence_name'] = seqnames
    summary['chromosome'] = chroms
    summary['start'] = starts
    summary['end'] = ends
    summary['sequence'] = seqs
    summary['score'] = scores
    summary['pvalue'] = pvalues
    summary['strand'] = strands
    
    returnDict[psid] = summary
    seenSeqsDict[psid] = seqsSeen
    seenNucsDict[psid] = seqsSeen * motif.getWidth() # in every sequence we have seen motif_width nucleotides


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


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
        Print the progress bar in the sequence scoring process
        ----
        Parameters:
            iteration (int)
            total (int)
            prefix (str)
            suffix (str)
            decimals (int)
            length (int)
            fill (str)
            printEnd (str)
        ----
        Returns
            None
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()   
        
def buildDF(motif, seqnames, starts, ends, strands, 
                scores, pvalues, qvalues, sequences, pvalueT):
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

    # all lists must have the same length
    LSTLEN=len(seqnames)

    assert len(starts) == LSTLEN
    assert len(ends) == LSTLEN
    assert len(strands) == LSTLEN
    assert len(scores) == LSTLEN
    assert len(pvalues) == LSTLEN
    assert len(sequences) == LSTLEN

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

            if QVAL:
                qvalues_thresh.append(qvalues[idx])

    DFLEN = len(seqnames_thresh)

    motifIDs=[motif.getMotifID()]*DFLEN
    motifNames=[motif.getMotifName()]*DFLEN

    # build the final dataframe
    df=pd.DataFrame()
    df['motif_id']=motifIDs
    df['motif_alt_id']=motifNames
    df['sequence_name']=seqnames_thresh
    df['start']=starts_thresh
    df['stop']=ends_thresh
    df['strand']=strands_thresh
    df['score']=scores_thresh
    df['p-value']=pvalues_thresh

    # ad the q-values to the final dataframe if they have been computed
    if QVAL:
        df['q-value']=qvalues_thresh

    df['matched_sequence']=sequences_thresh
        
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

    print()
    for _ in range(20):
            print('#', end='')
    print() # newline
    print('Scoring hits for motif', posID, sep=' ')
    print()
    
    # if we score also the reverse complement
    if not no_reverse: 
        print('Scoring hits for motif', revID, sep=' ')
        print()
        
    for _ in range(20):
        print('#', end='')
    print() # newline
        