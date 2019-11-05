"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Script that performs the extraction of the subgraphs from the graph genome and 
takes the sequences to score from them

"""

import subprocess
import warnings
from . import handle_exception as he
import sys
import os
import numpy as np
import multiprocessing as mp

peak_no=int(1) # global var

def get_data(genome_loc, bedfile, TFBS_len, vg_creation_pipeline, gplus, 
                 chroms, ncores):
    """
        returns the data necessary for the following steps of the analysis 
        pipeline
        ----
        Parameters:
            genome_loc (str) : path to the genome 
            bedfile (str) : path to the bedfile
            TFBS_len (int) : width of the motif 
            vg_creation_pipeline (bool) : defines which pipeline to follow
        ----
        Returns:
            fileloc (str) : path to the data obtained
    """
    
    for _ in range(20):
        print('#', end='')
    print() #newline
    print("Extracting subgraphs from regions defined in ", bedfile)
    print()
    for _ in range(20):
        print('#', end='')
    print()
    
    tmpwd='.grafimo'
    cmd='mkdir -p {0}'.format(tmpwd)
    code=subprocess.call(cmd, shell=True)
    if code!=0:
        raise Exception('error while executing a shell command (mkdir)')
        sys.exit(1)
    
    try:
        
        NCORES=ncores
        bedregions=getBEDregions(bedfile)
        
        # prepare for parellelization
        jobs=[]
        bedregions_splt=np.array_split(bedregions, NCORES)
        
        cwd=os.getcwd()
        os.chdir(tmpwd)
        
        if bedfile[0:6]=='/Users':
            pass # given an absolute path --> no problems
        else:
            bedfile=''.join(['../', bedfile]) # given a relative path
                                              # NB we are inside .grafimo
        
        if vg_creation_pipeline:
            
            for i in range(NCORES):
                p=mp.Process(target=vgc_sge, args=(bedregions_splt[i], 
                                                       TFBS_len, chroms))
                jobs.append(p)
                p.start()
                
            for job in jobs:
                job.join() # deadlock
        
        elif not vg_creation_pipeline:
    
            if gplus:
                for i in range(NCORES):
                    p=mp.Process(target=no_vgc_sge_gplus, 
                                     args=(genome_loc, bedregions_splt[i], 
                                               TFBS_len, chroms))
                    jobs.append(p)
                    p.start()
                    
                for job in jobs:
                    job.join() # deadlock
        
            else:
                if genome_loc[0]=='~':
                    pass
                else:
                    genome_loc=''.join(['../', genome_loc])
                    
                for i in range(NCORES):
                    p=mp.Process(target=no_vgc_sge, 
                                     args=(genome_loc, bedregions_splt[i],
                                               TFBS_len))
                    jobs.append(p)
                    p.start()
                    
                for job in jobs:
                    job.join() # deadlock
            
        else:
            # with enter an unknown value for the pipeline
            raise Exception('Error: pipeline not recognized')
            sys.exit(1) # the workframe fails
        
    except:
        cmd='rm -rf {0}'.format(tmpwd) # clean the directory from tmp data
        code=subprocess.call(cmd)
        if code!=0:
            raise Exception('Error: error while executing a shell command (rm)')
            sys.exit(1) # exit for the shell cmd error
       
        raise Exception('Error: subgraphs extraction failed')
        sys.exit(1) # failure in subgraphs extraction
        
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
        if toXGpath[-i]=='/':
            bp=-i
            break
        
    toXGpath=toXGpath[:bp]
    
    return toXGpath

def vgc_sge(bedfile, TFBS_len, chroms):
    """
        Extract the subgraphs from the genome graph and the sequences to score
        from them (step to follow in the pipeline with the vg creation)
        ----
        Parameters:
            bedfile (str) : path to the bedfile
            TFBS_len (int) : motif width
        ----
        Returns:
            None
    """
    
#    if bedfile[0:6]=='/Users':
#        pass
#    else:
#        bedfile=''.join(['../', bedfile])
    
    
    try:
        #bed=open(bedfile, mode='r') # open the bedfile in read only mode
        
        global peak_no
        
        if chroms:
            CHR_LIST=[''.join(['chr', c]) for c in chroms]
            
        else:
            CHR_LIST=list(range(1,23))+['X', 'Y']
            CHR_LIST=[''.join(['chr', str(c)]) for c in CHR_LIST]
        
        for line in bedfile:
            print('Peak ', peak_no)
            chrom, start, end = line.split('\t')[0:3]
            print('Extracting region:', chrom, '['+start+'-'+end+']')
            #chromid=chrom.replace('chr', '')
            
            if chrom in CHR_LIST: # chromosome name is valid
                region_index=chrom+':'+start+'-'+end
            
                path_id=chrom+'_'+start+'-'+end
                subgraph_path=correct_path('./', path_id, '.vg')
            
                vg=''.join(['../', chrom, '.xg']) # for extraction xg is required
                
                code=extract_region(vg, region_index, subgraph_path)
            
                if code != 0:
                    warn="Region "+chrom+':'+'['+start+'-'+end+']'+" extraction failed!\n"
                    warnings.warn(warn, Warning)  # although we have an exception we don't stop the execution
                    
                else:
                    # write on the stderr
                    msg="Region "+chrom+':['+start+'-'+end+']'+" extracted\n"
                    sys.stderr.write(msg) 
                
                kmer_path=correct_path('./', path_id, '.tsv')
                
                # vg kmer works with vg not xg graphs
                code=retrieve_kmers(TFBS_len, subgraph_path, kmer_path)
            
                if code != 0:
                    warn="Kmers extraction from region "+chrom+':'+'['+start+'-'+end+']'+" failed!\n"
                    warnings.warn(warn, Warning) # although we have an exception we don't stop the execution
                    
                else:
                    # write on the stderr
                    msg="Kmers extraction for region "+chrom+':['+start+'-'+end+']'+" finished\n"
                    sys.stderr.write(msg)
                
                cmd='rm {0}'.format(subgraph_path)
                code=subprocess.call(cmd, shell=True)
            
                if code != 0:
                    cderr=he.throw_subprocess_error()
                    sys.exit(cderr)
                    
            else:
                warnings.warn('chromosome name not valid. Region skipped\n', Warning)
                # although we have an exception we don't stop the execution
                          
            peak_no += 1
            
    except IOError:
        sys.stderr.write('Error: unable to open', bedfile) 
        # something went wrong while reading the bedfile
        sys.exit(1)
        
    #else:
     #   bed.close() # close the bedfile
        

def no_vgc_sge(xg, bedfile, TFBS_len):
    """
        Extract the subgraphs from the genome graph and the sequences to score
        from them (step to follow in the pipeline without the vg creation)
        ----
        Parameters:
            bedfile (str) : path to the bedfile
            TFBS_len (int) : motif width
        ----
        Returns:
            None
    """
#    if bedfile[0:6]=='/Users':
#        pass # given an absolute path
#    else:
#        bedfile=''.join(['../', bedfile])
    
    global peak_no
    
    try:
        #bed=open(bedfile, mode='r') # open the bedfile in read only mode
        
        CHR_LIST=list(range(1,23))+['X', 'Y']
        CHR_LIST=[''.join(['chr', str(c)]) for c in CHR_LIST]
        
        for line in bedfile:
            print('Peak ', peak_no)
            chrom, start, end = line.split('\t')[0:3]
            print('Extracting region:', chrom, '['+start+'-'+end+']')
            #chromid=chrom.replace('chr', '')
            
            if chrom in CHR_LIST: # chromosome name is valid
                region_index=chrom+':'+start+'-'+end
            
                path_id=chrom+'_'+start+'-'+end
                subgraph_path=correct_path('./', path_id, '.vg')
                
                code=extract_region(xg, region_index, subgraph_path)
            
                if code != 0:
                    warn="Region "+chrom+':'+'['+start+'-'+end+']'+" extraction failed!\n"
                    warnings.warn(warn, Warning)  # although we have an exception we don't stop the execution
                    
                else:
                    # write to stderr
                    msg="Region "+chrom+':['+start+'-'+end+']'+" extracted\n"
                    sys.stderr.write(msg) 
                
                kmer_path=correct_path('./', path_id, '.tsv')
                
                # vg kmer works with vg not xg graphs
                code=retrieve_kmers(TFBS_len, subgraph_path, kmer_path)
            
                if code != 0:
                    warn="Kmers extraction from region "+chrom+':'+'['+start+'-'+end+']'+" failed!\n"
                    warnings.warn(warn, Warning) # although we have an exception we don't stop the execution
                    
                else:
                    # write on the stderr
                    msg="Kmers extraction for region "+chrom+':['+start+'-'+end+'] '+"finished\n"
                    sys.stderr.write(msg)
                
                cmd='rm {0}'.format(subgraph_path)
                code=subprocess.call(cmd, shell=True)
            
                if code != 0:
                    cderr=he.throw_subprocess_error()
                    sys.exit(cderr)
                    
            else:
                warnings.warn('chromosome name not valid. Region skipped\n', Warning)
                # although we have an exception we don't stop the execution
                          
            peak_no += 1
            
    except IOError:
        print('Error: unable to open', bedfile) 
        # something went wrong while reading the bedfile
        sys.exit(1)
        
    #finally:
        #bed.close() # close the bedfile
        

def no_vgc_sge_gplus(xg, bedfile, TFBS_len, chroms):
        
    # check the path to xgs
    if xg[0:6]=='/Users':
        pass
    else:
        xg=''.join(['../', xg])
    
    global peak_no
    
    try:
        #bed=open(bedfile, mode='r')
        
        if chroms:
            CHR_LIST=[''.join(['chr', c]) for c in chroms]
            
        else:
            CHR_LIST=list(range(1,23))+['X', 'Y']
            CHR_LIST=[''.join(['chr', str(c)]) for c in CHR_LIST]
        
        for line in bedfile:
            
            print('Peak ', peak_no)
            chrom, start, end = line.split('\t')[0:3]
            print('Extracting region:', chrom, '['+start+'-'+end+']')
            
            if chrom in CHR_LIST: # chromosome name is valid
                region_index=chrom+':'+start+'-'+end
            
                path_id=chrom+'_'+start+'-'+end
                subgraph_path=correct_path('./', path_id, '.vg')
                
                vg=''.join([xg, chrom, '.xg']) # for extraction xg is required
                    
                code=extract_region(vg, region_index, subgraph_path)
            
                if code != 0:
                    warn="Region "+chrom+':'+'['+start+'-'+end+']'+" extraction failed!\n"
                    warnings.warn(warn, Warning)  # although we have an exception we don't stop the execution
                    
                else:
                    # write on the stderr
                    msg="Region "+chrom+':['+start+'-'+end+']'+" extracted\n"
                    sys.stderr.write(msg) 
                
                kmer_path=correct_path('./', path_id, '.tsv')
                    
                # vg kmer works with vg not xg graphs
                code=retrieve_kmers(TFBS_len, subgraph_path, kmer_path)
                    
                if code != 0:
                    warn="Kmers extraction from region "+chrom+':'+'['+start+'-'+end+']'+" failed!\n"
                    warnings.warn(warn, Warning) # although we have an exception we don't stop the execution
                    
                else:
                    # write on the stderr
                    msg="Kmers extraction for region "+chrom+':['+start+'-'+end+']'+" finished\n"
                    sys.stderr.write(msg)
                
                cmd='rm {0}'.format(subgraph_path)
                code=subprocess.call(cmd, shell=True)
            
                if code != 0:
                    cderr=he.throw_subprocess_error()
                    sys.exit(cderr)
                    
            else:
                warnings.warn('chromosome name not valid. Region skipped\n', Warning)
                # although we have an exception we don't stop the execution
                          
            peak_no += 1
            
    except IOError:
        print('Error: unable to open', bedfile) # something went wrong while reading the bedfile
        sys.exit(1)
        
    else:
        #bed.close() # close the bedfile
        pass
        
    
def correct_path(path, path_id='', file_format=''):
    
    if path[-1:]=='/':
        new_path=path+path_id+file_format
    else:
        new_path=path+'/'+path_id+file_format
        
    return new_path

def retrieve_kmers(TFBS_len, subgraphs_path, kmers_path):
    """
        Obtain the sequences from the subgraphs
        ----
        Prameters:
            TFBS_len (int) : motif width
            subgraphs_path (str) : path to the subgraphs
            kmers_path (str) : path to which the sequences will be stored
        ----
        Returns:
            code (int) : success or not of subprocess.call()
    """    
    vg_km_cmd='vg kmers -k {0} -p {1} > {2}'.format(str(TFBS_len), subgraphs_path, kmers_path)
    code=subprocess.call(vg_km_cmd, shell=True)
    
    return code

def extract_region(genome, region_index, subgraphs_path):
    """
        Extract the subgraphs
        ----
        Parameters:
            genome (str) : path to the genome
            region_index (str) : region to extract
            subgraphs_path (str) : path to which the subgraphs will be stored
        ----
        Returns:
            code (int) : success or not of subprocess.call()
    """
    
    vg_sg_cmd='vg find -x {0} -p {1} > {2}'.format(genome, region_index, subgraphs_path)
    code=subprocess.call(vg_sg_cmd, shell=True)
    
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
        code=he.throw_not_str_error()
        sys.exit(code)
        
    else:
        if graph_genome.split()[-1] == 'xg':
            return True
        else:
            return False
    
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
        raise Exception("Error: the given file is not in BED format")
        sys.exit(1)
        
    regions=[]    
    
    try:
        b=open(bedfile, mode='r') # open the BED file in read only mode
        
        for line in b:
            regions.append(line)
            
    except: # not able to read the BED file
        msg=' '.join(["Error: unable to read", bedfile])
        raise Exception(msg)
        sys.exit(1)
        
    else:
        return regions
    
    