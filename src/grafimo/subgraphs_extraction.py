"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Script that performs the extraction of the subgraphs from the graph genome and 
takes the sequences to score from them

"""

import subprocess
import warnings
from grafimo.GRAFIMOException import SubprocessException, NotValidFFException, FileReadingException, VGException, \
                                        PipelineException, ValueException
from grafimo.utils import die, correct_path, CHROMS_LIST, printProgressBar
import sys
import os
import numpy as np
import multiprocessing as mp
import tempfile

####################################################
#
#TODO: support also non conventional chromosomes
#
####################################################

def get_data(genome_loc, bedfile, TFBS_len, vg_creation_pipeline, gplus, 
                 chroms, ncores, verbose = False):
    """
        returns the data necessary for the following steps of the analysis 
        pipeline
        ----
        Parameters:
            genome_loc (str) : path to the genome 
            bedfile (str) : path to the bedfile
            TFBS_len (int) : width of the motif 
            vg_creation_pipeline (bool) : defines which pipeline to follow
            gplus (bool) : flag value that says if we have different chromosome
                           graphs or a single genome graph
            chroms (list) : chromosomes from which regions will be extracted
            ncores (int) : number of cores to use
        ----
        Returns:
            fileloc (str) : path to the data obtained
    """

    printSgeWelcomeMsg(bedfile)

    # create a tmp working directory
    tmpwd = tempfile.mkdtemp(prefix = 'grafimo_')

    # if the tmp directory name already exists remove it
    # this shouldn't happen, but to be sure...
    if os.path.isdir(tmpwd):
        cmd='rm -rf {0}'.format(tmpwd)
        code=subprocess.call(cmd, shell=True)

        if code != 0:
            raise SubprocessException(' '.join(["an error occurred executing", cmd, ". Exiting"]))
            die(1)

    cmd = 'mkdir -p {0}'.format(tmpwd)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise SubprocessException(' '.join(["an error occurred executing", cmd, ". Exiting"]))
        die(1)
    
    try:
        
        NCORES = ncores
        bedregions = getBEDregions(bedfile)
        
        # prepare BED for parellelization
        jobs = []
        proc_finished = 0
        bedregions_splt = np.array_split(bedregions, NCORES)

        # get the new location of graphs wrt the tmp dir
        cwd = os.getcwd()
        os.chdir(genome_loc)
        genome_loc = os.getcwd() # now we are sure that it's an absolute path
        os.chdir(cwd) # get back to the starting point 

        # enter the tmp dir
        os.chdir(tmpwd)

        # print progress bar if the verbose == False
        if not verbose:
            printProgressBar(proc_finished, NCORES, prefix='Progress:',
                                    suffix='Complete', length=100)

        if vg_creation_pipeline:
            
            for i in range(NCORES):
                p=mp.Process(target=vgc_sge, args=(genome_loc, bedregions_splt[i], 
                                                       TFBS_len, chroms, verbose,cwd))
                jobs.append(p)
                p.start()
                
            for job in jobs:
                proc_finished += 1
                
                if not verbose:
                    printProgressBar(proc_finished, NCORES, prefix='Progress:',
                                        suffix='Complete', length=100)
                
                job.join() # deadlock
        
        elif not vg_creation_pipeline:
    
            if gplus: # we have more chromosome graphs
                for i in range(NCORES):
                    p=mp.Process(target=no_vgc_sge_gplus, 
                                     args=(genome_loc, bedregions_splt[i], 
                                               TFBS_len, chroms, verbose, cwd))
                    jobs.append(p)
                    p.start()
                    
                for job in jobs:
                    proc_finished += 1
                    
                    if not verbose:
                        printProgressBar(proc_finished, NCORES, prefix='Progress:',
                                            suffix='Complete', length=100)
                    
                    job.join() # deadlock
        
            else: # we have a single genome graph
                    
                for i in range(NCORES):
                    p=mp.Process(target=no_vgc_sge, 
                                     args=(genome_loc, bedregions_splt[i],
                                               TFBS_len, verbose, cwd))
                    jobs.append(p)
                    p.start()
                    
                for job in jobs:
                    proc_finished += 1
                    
                    if not verbose:
                        printProgressBar(proc_finished, NCORES, prefix='Progress:',
                                            suffix='Complete', length=100)
                    
                    job.join() # deadlock
            
        else:
            # with enter an unknown value for the pipeline
            raise PipelineException('Pipeline not recognized. Cannot proceed')
            die(1) # this should never happen
        
    except:
        cmd='rm -rf {0}'.format(tmpwd) # clean the directory from tmp data
        code=subprocess.call(cmd)

        if code!=0:
            raise SubprocessException(' '.join(["an error occurred executing", cmd, ". Exiting"]))
            die(1) # exit for the shell cmd error
       
        raise VGException('Subgraphs extraction failed')
        die(1) # failure in subgraphs extraction
        
    else:
        
        fileloc=os.getcwd() # the sequences to score are store in the cwd 
    
        os.chdir(cwd) # get back to the original directory
        
        return fileloc

def get_xg_loc(toXGpath):
    """
        Get the path to the directory that contains the whole genome XG
        ----
        Parameters:
            toXGpath(str) : path to the xg
        ----
        Returns:
            toXGpath (str) : path to the directory containing the xg
    """
    
    for i in range(1, len(toXGpath)):
        if toXGpath[-i] == '/':
            bp = -i
            break
        
    toXGpath = toXGpath[:bp]
    
    return toXGpath

def vgc_sge(xg_loc, bedfile, TFBS_len, chroms, verbose, cwd):
    """
        Extract the subgraphs from the genome graph and the sequences to score
        from them (step to follow in the pipeline with the vg creation)
        ----
        Parameters:
            xg_loc (str) : path to the genome graphs location
            bedfile (list) : slice of regions to extract
            TFBS_len (int) : motif width
            chroms (list) : list of chromosomes to take into account during the
                            extraction
        ----
        Returns:
            None
    """
    
    if xg_loc[-1] == "/":
        pass
    else:
        xg_loc = ''.join([xg_loc, "/"])
    
    try:

        # look at user defined chromosomes
        if chroms:
            CHR_LIST = [''.join(['chr', c]) for c in chroms]

        # look at all chromosomes
        else:
            CHR_LIST = [''.join(['chr', c]) for c in CHROMS_LIST]
        
        for line in bedfile:
            chrom, start, end = line.split('\t')[0:3]
            #chromid=chrom.replace('chr', '')
            
            if chrom in CHR_LIST: # chromosome name is valid
                region_index = chrom + ':' + start + '-' + end
                if verbose:
                    print('Extracting region:', region_index)
            
                path_id = chrom + '_' + start + '-' + end
                kmers_file = correct_path('./', path_id, '.tsv')
            
                xg = ''.join([xg_loc, chrom, '.xg']) # for extraction is required the xg
                
                code = get_kmers(xg, region_index, TFBS_len, kmers_file)
            
                if code != 0:
                    if verbose:
                        warn = "Region "+chrom+':'+'['+start+'-'+end+']'+" extraction failed!\n"
                        sys.stderr.write(warn)  # it can happen that wasn't possible to extract a peak
                                                # but we don't stop the execution
                    else:
                        pass
                    
                else:
                   
                    if verbose:
                        # write on the stderr
                        msg = "Region "+chrom+':['+start+'-'+end+']'+" extracted\n"
                        sys.stderr.write(msg)
                    else:
                        pass 
                
                # if was not possible to extract the region we go to the next one
                if os.stat(kmers_file).st_size <= 0:

                    cmd = 'rm {0}'.format(kmers_file) # remove the empty file 
                    code = subprocess.call(cmd, shell=True)
            
                    if code != 0:
                        raise SubprocessException(' '.join(["An error occurred while executing", cmd, ". Exiting"]))
                        die(1)
                    
                    continue
                    
            else:
                if verbose:
                    warn = 'chromosome name not valid. Region skipped\n'
                    sys.stderr.write(warn) # only canonical chromosomes are supported
                else:   
                    pass
            
    except:
        raise FileReadingException("Unable to extract regions")
        die(1)

    else:
        pass
        

def no_vgc_sge(xg, bedfile, TFBS_len, chroms, verbose, cwd):
    """
        Extract the subgraphs from the genome graph and the sequences to score
        from them (step to follow in the pipeline without the vg creation)
        ----
        Parameters:
            xg (str) : genome graph
            bedfile (str) : path to the bedfile
            TFBS_len (int) : motif width
            chroms (list) : list of chromosomes from which the regions
                            will be extracted
        ----
        Returns:
            None
    """
    
    try:

        if chroms:
            CHR_LIST = [''.join(['chr', c]) for c in chroms]
            
        else:
            CHR_LIST = [''.join(['chr', c]) for c in CHROMS_LIST]
        
        for line in bedfile:

            chrom, start, end = line.split('\t')[0:3]
            
            if chrom in CHR_LIST: # chromosome name is valid
                region_index = chrom + ':' + start + '-' + end
                
                if verbose:
                    print("Extracting region:", region_index)
            
                path_id = chrom + '_' + start + '-' + end
                kmers_file = correct_path('./', path_id, '.tsv')
                
                code = get_kmers(xg, region_index, TFBS_len, kmers_file)
            
                if code != 0:
                    if verbose:
                        warn = "Region " + chrom + ':' + '[' + start + '-' + end + ']' + " extraction failed!\n"
                        sys.stderr.write(warn)  # it can happen that wasn't possible to extract a peak
                                                # but we don't stop the execution
                    else:
                        pass

                else:
                    if verbose:
                        # write to stderr
                        msg = "Region " + chrom + ':[' + start + '-' + end + ']' + " extracted\n"
                        sys.stderr.write(msg)
                    else:
                        pass 

                # if was not possible to extract the region we go to the next one
                if os.stat(kmers_file).st_size <= 0:

                    cmd = 'rm {0}'.format(kmers_file) # remove the empty file 
                    code = subprocess.call(cmd, shell=True)
            
                    if code != 0:
                        raise SubprocessException(' '.join(["An error occurred while executing", cmd, ". Exiting"]))
                        die(1)
                    
                    continue
                    
            else:
                if verbose:
                    warn = 'chromosome name not valid. Region skipped\n'
                    sys.stderr.write(warn) # only canonical chromosomes are supported
                else:
                    pass
            
    except:
        raise FileReadingException("Unable to extract regions")
        die(1)
        
    else:
        pass

def no_vgc_sge_gplus(xg_loc, bedfile, TFBS_len, chroms, verbose, cwd):
    """
        Extract the subgraphs from the selected chromosomes and
        the sequences to score from them (step to follow in the
        pipeline without the vg creation)
        ----
        Parameters:
            xg_loc (str) : path to the chromosomes variation graphs
            bedfile (list) : slice of regions to extract
            TFBS_len (int) : motif width
            chroms (list) : list of chromosomes from which the regions
                            will be extracted
        ----
        Returns:
            None
    """
    
    if xg_loc[-1] == '/':
        pass
    else:
        xg_loc = ''.join([xg_loc, '/'])

    try:
        
        if chroms:
            CHR_LIST = [''.join(['chr', c]) for c in chroms]
            
        else:
            CHR_LIST = [''.join(['chr', c]) for c in CHROMS_LIST]
        
        for line in bedfile:

            chrom, start, end = line.split('\t')[0:3]
            
            if chrom in CHR_LIST: # chromosome name is valid
                region_index = chrom + ':' + start + '-' + end
                
                if verbose:
                    print('Extracting region:', region_index)
            
                path_id = chrom + '_' + start + '-' + end
                kmers_file = correct_path('./', path_id, '.tsv')
                
                xg = ''.join([xg_loc, chrom, '.xg']) # for extraction is required the xg
                    
                code = get_kmers(xg, region_index, TFBS_len, kmers_file)
            
                if code != 0:
                    if verbose:
                        warn = "Region " + chrom + ':' + '[' + start + '-' + end + ']' + " extraction failed!\n"
                        sys.stderr.write(warn)  # it can happen that wasn't possible to extract a peak
                                                # but we don't stop the execution
                    else:
                        pass
                    
                else:
                    if verbose:
                        # write on the stderr
                        msg="Region "+chrom+':['+start+'-'+end+']'+" extracted\n"
                        sys.stderr.write(msg)
                    else:
                        pass 

                # if was not possible to extract the region we go to the next one
                if os.stat(kmers_file).st_size <= 0:

                    cmd = 'rm {0}'.format(kmers_file) # remove the empty file 
                    code = subprocess.call(cmd, shell=True)
            
                    if code != 0:
                        raise SubprocessException(' '.join(["An error occurred while executing", cmd, ". Exiting"]))
                        die(1)
                    
                    continue
                
            else:
                if verbose:
                    warn = 'chromosome name not valid. Region skipped\n'
                    sys.stderr.write(warn) # only canonical chromosomes are supported
                else:
                    pass # we simply ignore this situation
            
    except:
        raise FileReadingException("Unable to extract regions")
        die(1)
        
    else:
        pass


def get_kmers(xg, region, TFBS_len, seqs_file):
    """
        Obtain the sequences of length k inside the queried
        region
        ----
        Prameters:
            xg (str) : path to the 
            region (str) : region to query on the graph
            TFBS_len (int) : motif width
            seqs_file (str) : name of the file where the
                                retrieved k-mers will be
                                put
        ----
        Returns:
            code (int) : success or not of subprocess.call()
    """    
    cmd = 'vg find -x {0} -E -p {1} -K {2} > {3}'.format(xg, region, TFBS_len, seqs_file)
    code = subprocess.call(cmd, shell=True)
    
    return code


def isGraph_genome_xg(graph_genome):
    """
        Check if the given genome graph is in xg format
        ----
        Parameters:
            graph_genome (str) : path to the genome graph
        ----
        Returns:
            (bool)
    """
    
    if not isinstance(graph_genome, str):
        raise VGException("Invalid path to the genome graph. Cannot proceed")
        die(1)

    if graph_genome.split('.')[-1] == 'xg':
        return True
    elif graph_genome.split('.') == 'vg':
        return False
    else:
        msg="Do not know what to do with the given genome graph. Only XG or VG format allowed"
        raise VGException(msg)
        die(1)
    

def getBEDregions(bedfile):
    """
        Function to retrieve the number of regions defined in the given BED
        ----
        Parameters:
            bedfile (str) : path to the bedfile to read
        ----
        Returns: 
            regions (list) : array of the regions defined in the BED file
    """
    
    if bedfile.split('.')[-1]!='bed': #not a BED file
        raise NotValidFFException("The given BED file is not in BED format")
        die(1)
        
    regions=[]    
    
    try:
        b=open(bedfile, mode='r') # open the BED file in read only mode
        
        for line in b:
            regions.append(line)
            
    except: # not able to read the BED file
        msg=' '.join(["Error: unable to read", bedfile])
        raise FileReadingException(msg)
        die(1)
        
    else:
        return regions

    finally:
        b.close() # close the file stream

def printSgeWelcomeMsg(bedfile):

    print()  # newline
    print("\nExtracting subgraphs from regions defined in ", bedfile)
    print()

