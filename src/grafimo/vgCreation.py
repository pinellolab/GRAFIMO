"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Construction and indexing of the genome-graphs.
The reference used is given by the user as like as the VCF 
used to infer the variants and the alternative paths on 
the graphs.

We choose to split the genome graph for the genome in a graph for
each chromosome, in order to reduce the amount of memory and
computational resources needed.

NB There are NO drawbacks using this approach instead of creating the whole genome graph

"""

import subprocess
import os
import glob
from grafimo.utils import CHROMS_LIST, die
from grafimo.GRAFIMOException import VGException, ValueException, SubprocessException


def get_genome_from_ucsc():
    """
        Download the reference genome hg38 from UCSC database and
        returns the path to it
        ----
        Parameters:
            None
        ----
        Returns:
            genome (str) : path to the genome downloaded (in .fa format)
    
    """
    
    # download genome
    cmd = 'wget -c ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
    code = subprocess.call(cmd, shell=True)
    
    if code != 0:
        raise SubprocessException(' '.join(["An error occurred while executing", cmd, ". Exiting"]))
        die(1)
        
    else:
        # decompress genome
        genome_comp = './hg38.fa.gz'
        genome = './hg38.fa'
        
        cmd = 'gunzip {0}'.format(genome_comp)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            raise SubprocessException(' '.join(["An error occurred while executing", cmd, ". Exiting"]))
            die(1)

        # remove .fasta.gz file
        cmd = 'rm {0}'.format(genome_comp)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            raise SubprocessException(' '.join(["An error occurred while executing", cmd, ". Exiting"]))
            die(1)
            
        else:
            return genome
        
    return os.getcwd() # if something goes wrong we return the cwd (this line doesn't convince me)


def get_1000GProject_vcf():
    """
        Downloads the vcf file from The 1000 Genome Project database and
        returns the path to it
        ----
        Parameters:
            None
        ----
        Returns:
            vcf (str) : path to the vcf downloaded vcf file (in .vcf.gz)
    """
    
    # download the VCF
    cmd = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/'
    cmd += '1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/'
    cmd += 'ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz'
    
    code = subprocess.call(cmd, shell=True)
    
    if code != 0:
        raise SubprocessException(' '.join(["An error occurred while executing", cmd, ". Exiting"]))
        die(1)
        
    else:
        vcf = './ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz'
        return vcf
            
        
def create_vg(chroms, linear_genome='', vcf=''):
    """
        Create the genome graph, for the reference and VCF file given
        in input by the user.

        The genome is not built as a single whole genome graph but a 
        single graph is constructed for each chromosome.
        This choice was made to avoid memory issues and make able
        also the less powerful machines to run GRAFIMO.

        There is NO drawback using this approach wrt  
        construct the whole genome graph and query it. 
        ----
        Parameters:
            chroms (list) : list of chromosomes for whicgh the genome
                            graph will be constructed
            linear_genome (str) : path to the linear genome used as
                                    reference to build the genome
                                    graphs
            vcf (str) : path to the VCF file used to build the genome
                        graphs
        ----
        Returns:
            graphs_loc (str) : path to the location of the genome graphs        
    """
    
    if not linear_genome:
        # get hg38 from UCSC
        linear_genome = get_genome_from_ucsc()
    elif linear_genome:
        NO_LINEAR_GENOME_INPUT = False
        
    if not vcf:
        # get VCF from The 1000 Genome Project
        vcf = get_1000GProject_vcf()
    elif vcf:
        NO_VCF_INPUT = False
        
    #if not chroms:
    #    chroms = CHROMS_LIST
        
    cwd = os.getcwd()
    
    if not tbiexist(cwd):
        msg = ''.join([vcf.split('/')[-1], ".tbi file not found indexing with tabix..."])
        print(msg)
        cmd = 'tabix -p vcf {0}'.format(vcf)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            raise SubprocessException(' '.join(["An error occurred while executing", cmd, ". Exiting"]))
            die(1)
    else: # we update the index
        msg = ' '.join(["Reindexing", vcf.split('/')[-1], "..."])
        print(msg)

        # remove the existing TBI file
        cmd = "rm {0}".format(''.join([vcf, ".tbi"]))
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            raise SubprocessException(' '.join(["An error occurred while executing", cmd, ". Exiting"]))
            die(1)

        # reindex the VCF
        cmd = "tabix -p vcf {0}".format(vcf)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            raise SubprocessException(' '.join(["An error occurred while executing", cmd, ". Exiting"]))
            die(1)
    
    for chrom_n in chroms:
        chrom = ''.join(['chr', chrom_n])
            
        vg = chrom+'.vg'
        xg = chrom+'.xg'
        
        """
        vg_construct = 'vg construct -C -R {0} -p -n {1}={2} -r {3} -v {4} > {5}'.format(chrom_n, chrom_n, \
                                                                                            chrom, linear_genome, vcf, vg)
        
        code=subprocess.call(vg_construct, shell=True)

        if code!=0:
            msg='error in vg construct. '
            msg+='Unable to build the vg of the genome using {0} and {1}'.format(linear_genome,
                                                                                    vcf)
            raise VGException(msg)
        """

        constructVG(vg, linear_genome, vcf, chrom, chrom_n)

        msg = ' '.join(["Indexing", vg, '...'])
        print(msg)
            
        """
        vg_index = 'vg index -p -x {0} {1}'.format(xg, vg)
        code = subprocess.call(vg_index, shell=True)

        if code != 0:
            msg = 'error in vg index. '
            msg += 'Unable to build the vg of the genome using {0} and {1}'.format(linear_genome,
                                                                               vcf)
            raise VGException(msg)
        """
        indexVG(vg)

        cmd = 'rm {0}'.format(vg)
        subprocess.call(cmd, shell=True)
        
        if code != 0:     # we have errors in the vg indexing
            raise SubprocessException(' '.join(["An error occurred while executing", cmd, ". Exiting"]))
            die(1)
            
    graphs_loc = os.getcwd()
    
    return graphs_loc

def constructVG(vg, linear_genome, vcf, chrom_name, chrom_num):
    """
        Build the genome graph of the given chromosome
        ----
        Parameters:
            vg (str) : output name of the genome graph
            linear_genome (str) : path to the linear genome location
            vcf (str) : path to the VCF file
            chrom_name (str) : chromosome name as used by UCSC (e.g. chr1)
            chrom_num (int) : chromosome number
        ----
        Returns:
            None
    """

    vg_construct = 'vg construct -C -R {0} -p -n {1}={2} -r {3} -v {4} > {5}'.format(chrom_num, chrom_num, \
                                                                                            chrom_name, linear_genome, vcf, vg)
        
    code = subprocess.call(vg_construct, shell=True)

    if code != 0:
        msg = 'error in vg construct. '
        msg += 'Unable to build the vg of the genome using {0} and {1}'.format(linear_genome,
                                                                                    vcf)
        raise VGException(msg)


def indexVG(vg):
    """
        Index the vg genome graph [vg => xg]. This step is required to 
        query regions on the graph in an efficient way.
        ----
        Parameters:
            vg (str) : path to the VG genome graph
        ----
        Returns:
            None
    """
    
    if not isinstance(vg, str):
        raise ValueException("Not valid path to the genome graph. Cannot proceed")
        
    else:
        xg = vg.split('.')[-2] # please, avoid paths that contains '.' in the middle
        xg += '.xg'
        vg_index = 'vg index -p -x {0} {1}'.format(xg, vg) 
        code = subprocess.call(vg_index, shell=True)
        
        if code != 0:
            raise VGException('Error in vg index. Unable to index the genome graph')
            die(1)
        

def tbiexist(directory):
    """
        Check for the existance of some indexed VCF file 
        (TBI file extension) in the specified directory
        ----
        Parameters:
            directory (str) : path to a directory
        ----
        Returns:
            (bool) 
    """

    path = ''.join([directory, '/*.tbi'])
    tbi = glob.glob(path)

    if len(tbi) > 0:
        return True
    
    return False
        
        