"""Functions to build and index the genome variation graph which will be
scanned for occurrences of the given motif(s).

The VG genome variation graph is built starting from a genome reference
then enriched with the given variants. 
After the its construction, the genome graph is indexed obtaining the 
indexed graph (XG) and a GBWT index, telling which are the path to 
follow within the graph to retrieve the genomes of the samples from
which the variants come from. So, each path in the genome correspond 
to a sample genome. The GBWT file is mandatory to recover this 
information.

The construction of the genome graph is made chromosome by chromosome,
in order to speed-up the construction and indexing operation and allow
an efficient parallel approach while scanning them for occurrences of
the searched motifs.

This approach is equivalent to construct a single VG, but allows to save 
space on disk and save memory when scanning it.
"""


from grafimo.utils import die, sigint_handler
from grafimo.workflow import BuildVG
from grafimo.GRAFIMOException import VGException, ValueException, \
    SubprocessError, FileReadingException
from typing import List
import subprocess
import signal
import time
import os


def get_reference_genome_from_ucsc() -> str:
    """Download the reference genome (hg38 assembly), from the UCSC
    database, in the current working directory and returns the path to 
    the corresponding FASTA file.

    This function has been written only for test purposes
    
    Parameters
    ----------

    Returns
    -------
    str
        path to the downloaded FASTA file (in .fa format)
    """

    cmd: str
    code: int
    errmsg: str

    # download genome
    address = 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
    cmd = 'wget -c {0}'.format(address)
    # the genome will be downloaded in the current directory
    code = subprocess.call(cmd, shell=True)  
    if code != 0:
        errmsg = ''.join(["\n\nERROR: an error occurred while executing ", 
                               cmd, ". Exiting"])
        raise SubprocessError(errmsg)

    # decompress genome
    print("Uncompressing the genome...")
    
    genome_comp: str = './hg38.fa.gz'
    if not os.path.exists(genome_comp):
        errmsg = ''.join(["\n\nERROR: ", genome_comp, " not found"])
        raise FileNotFoundError(errmsg)

    cmd = 'gunzip {0}'.format(genome_comp)
    code = subprocess.call(cmd, shell=True)

    if code != 0:
        errmsg = ''.join(["\n\nERROR: an error occurred while executing ", cmd, 
                          ". Exiting"])
        raise SubprocessError(errmsg)

    # remove FASTA.GZ file if still present
    if os.path.exists(genome_comp):
        cmd = 'rm {0}'.format(genome_comp)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            errmsg = ''.join(["\n\nERROR: an error occurred while executing ", 
                                   cmd, ". Exiting"])
            raise SubprocessError(errmsg)

    # get the path to the genome file
    genome_uncomp: str = "./hg38.fa"
    assert os.path.exists(genome_uncomp)
    genome: str = os.path.abspath(genome_uncomp)

    return genome

# end of get_reference_genome_from_ucsc()


def get_1000GProject_vcf() -> str:
    """Downloads a WGS VCF file from the 1000 Genome Project database
    (phase 3), containing SNVs and indels. The present file is used for 
    VG construction and graph indexing test purposes. 
    
    Since the variants present in this file are not phased, it cannot be 
    used to build the GBWT index and the corresponding haplotypes cannot
    be used. To use this features we must phase the VCF.

    Parameters
    ----------
    
    Returns
    -------
    str
        path to the downloaded VCF file (compressed)
    """

    address: str
    cmd: str
    code: int
    errmsg: str

    # download the VCF
    address = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/'
    address += '1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/'
    address += 'ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz'

    cmd = 'wget -c {0}'.format(address)

    code = subprocess.call(cmd, shell=True)

    if code != 0:
        errmsg = ''.join(["\n\nERROR: An error occurred while executing ", 
                               cmd, ". Exiting"])
        raise SubprocessError(errmsg)

    vcf_file: str = './ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz'
    vcf: str = os.path.abspath(vcf_file)

    return vcf

# end of get1000GProject_vcf()


def construct_vg(buildvg_args: BuildVG) -> None:
    """ Create the genome graph from the given genome reference and 
    phased VCF file given.
    
    The genome is not built as a single whole genome graph but a
    single graph is constructed for each chromosome.
    This approach avoids memory issues and allows the genome variation
    graph construction also on machines with less resources.

    There is NO drawback using this approach with respect to build
    a whole genome graph and query it.
    
    Moreover, it allows parallel queries on the different chromosomes to
    be perfromed also on regular laptops (>= 16 GB of memory), which is 
    very difficult with a whole genome graph, that requires the user
    to set appropriately the number of cores to use. Anyway a whole
    genome graph can be queried using a regular laptop using one core.
    
    Parameters
    ----------
        buildvg_args : BuildVG
            container for the arguments required to build the genome 
            variation graph 
    """
    errmsg: str
    if not isinstance(buildvg_args, BuildVG):
        errmsg = "Unknown arguments object type. "
        errmsg += "Cannot build the genome variation graph. Exiting"
        raise ValueError(errmsg)
    
    # read the arguments to build the VGs
    reindex: bool = buildvg_args.get_reindex()
    chroms: List[str] = buildvg_args.get_chroms()
    threads: int = buildvg_args.get_cores()
    outdir: str = buildvg_args.get_outdir()
    verbose: bool = buildvg_args.get_verbose()
    test: bool = buildvg_args.get_test()  # manually set in the code

    reference: str
    vcf: str

    if test:
        reference = get_reference_genome_from_ucsc()
        vcf = get_1000GProject_vcf()

    else:
        reference = buildvg_args.get_reference_genome()
        vcf = buildvg_args.get_vcf()
    # end if

    if verbose:
        print("using reference genome: ", reference)
        print("Using VCF file: ", vcf, "\n\n")
    # end if

    msg: str

    if verbose:
        start_c: float = time.time()
        msg = "Reading chromosome from reference file for which will be built the VG..."
        print(msg)
        
    # read the chromosome present in the used reference file
    chroms_available: List[str] = get_chromlist(reference)
        
    if verbose:
        end_c: int = time.time()
        print("done in %.2fs" % (end_c - start_c))
        print("Found chromosomes:\n", chroms_available, end="\n\n")

    if len(chroms) == 1 and chroms[0] == 'ALL_CHROMS':
        chroms: List[str] = chroms_available
    else:
        # check if the given chromosome number is among those whose sequence 
        # is in the reference file
        if (any(True for c in chroms if c not in chroms_available)):
            raise ValueError("Unknown chromosome given")
        
    # end if

    cwd: str = os.getcwd()

    cmd: str
    code: int

    # check if the VCF file has already been indexed with tabix
    if not tbiexist(vcf):
        msg = ''.join(["TBI file not found for ", vcf.split('/')[-1], 
                       ". Indexing the VCF file with tabix..."])
        print(msg)
        cmd = 'tabix -p vcf {0}'.format(vcf)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            # tabix didn't work
            errmsg = ''.join(["\n\nERROR: an error occurred while executing ", 
                              cmd, ". Exiting"])
            raise SubprocessError(errmsg)

    elif reindex:  # the user want to reindex the VCF file with tabix
        msg = ''.join(["Reindexing ", vcf.split('/')[-1], "...\n"])
        print(msg)

        # remove the existing TBI file
        cmd = "rm {0}".format(''.join([vcf, ".tbi"]))
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            errmsg = ''.join(["\n\nERROR: an error occurred while executing ", 
                              cmd, ". Exiting"])
            raise SubprocessError(errmsg)

        # reindex the VCF
        cmd = "tabix -p vcf {0}".format(vcf)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            # tabix didn't work
            errmsg = ''.join(["\n\nERROR: an error occurred while executing ", 
                              cmd, ". Exiting"])
            raise SubprocessError(errmsg)
        # end if
    # end if

    # enter the output directory
    os.chdir(outdir)

    # build the VG for each chromosome or a user defined
    # subset of them
    for chrom_n in chroms:

        chrom: str = ''.join(['chr', chrom_n])
        vg: str = ''.join([".", chrom, '.vg'])

        # build the VG for the current chromosome
        if verbose:
            start_build: float = time.time()

        code = build_vg(vg, reference, vcf, chrom_n, threads)
        if code != 0:
            msg = '\n\nERROR: an error occurred during {0} construction. '.format(vg)
            msg += 'Unable to build the VG of the genome using {0} and {1}'.format(reference,
                                                                                   vcf)
            raise VGException(msg)
        # end if

        if verbose:
            end_build: float = time.time()
            msg = "Elapsed time to build {0} ".format(vg)   
            print(msg, "%.2fs" % (end_build - start_build))
        # end if

        # to query efficiently the VGs we index them (VG -> XG)
        if verbose:
            start_index: float = time.time()

        msg = ''.join(["Indexing ", vg, ' and building the GBWT index...'])
        print(msg)

        code = indexVG(vg, vcf, threads, verbose)

        if code != 0:
            errmsg = "\n\nERROR: an error occurred while indexing {0}.".format(vg)
            errmsg += "\nUnable to index {0}. Exiting".format(vg)
            raise VGException(errmsg)
        # end if

        if verbose:
            end_index: float = time.time()
            msg = "Elapsed time to index {0}".format(vg)
            print(msg, "%.2fs" % (end_index - start_index))
        # end if

        # The majority of applications work only with indexed graph,
        # so to save disk space is worth to delete the VGs and keep
        # only the XGs (is simple to get back using VG built-in functions)
        if verbose:
            print("Deleting {0}".format(vg))

        cmd = 'rm {0}'.format(vg)
        subprocess.call(cmd, shell=True)

        if code != 0:  # we have errors in the vg indexing
            errmsg = ''.join(["\n\nERROR: an error occurred while executing ", cmd, ". Exiting"])
            raise SubprocessError()
        # end if
    # end for

    # get the VGs location
    graphs_loc: str = os.getcwd()

    # return to the original working directory
    os.chdir(cwd)

# end of construct_vg()


def build_vg(
    vg: str, 
    ref_genome: str, 
    vcf: str, 
    chrom_num: int, 
    threads: int
) -> int:
    """Build the genome variation graph for the given chromosome. The 
    chromosome sequence will be enriched with the genomic variants
    present in the input phased VCF file.
        
    NB. to build the VG we need that the chromosome names on the FASTA
    file of reference genome have format '>1' and not '>chr1'. Otherwise
    when indexing it the haplotype linking cannot be made.

    TODO: fix chromosome naming linking when building GBWT
    
    Parameters
    ----------
    vg : str
        name of the resulting genome variation graph
    ref_genome : str
        path to the reference genome FASTA file
    vcf : str
        path to the phased VCF file (compressed)
    chrom_num : int
        chromosome number for which the genome variation graph will be 
        constructed (e.g. chromosome 1 ==> 1)
    threads : int 
        number of threads to use while constructing the current VG
        
    Returns
    -------
    int  
        status of VG construction (0 = all ok; 1 = an error occurred)
    """

    success: int
    build: str

    build = 'vg construct -t {4} -r {1} -v {2} -R {0} -C -a -p > {3}'.format(chrom_num,
                                                                         ref_genome, vcf,
                                                                         vg, threads)

    code: int = subprocess.call(build, shell=True)

    if code != 0:
        # an error occurred during the construction of the VG for
        # the current chromosome
        success = 1
    else:
        # all went well
        success = 0

    return success

# end of build_vg()


def indexVG(
    vg: str, 
    vcf: str, 
    threads: int, 
    verbose: bool
) -> int:
    """Construct the XG and GBWT indexes for the given genome variation 
    graph. These indexes are required to query the genome when extracting
    motif occurrence candidates.

    The GBWT index allows to keep track of the haplotypes used to build
    the graph data structure and retrieve the samples genomes. 
    
    The indexing operation could take some time.
        
    Parameters
    ----------
    vg : str
        path to the genome variation graph (VG format)
    vcf : str 
        path to the phased VCF file used to build the corresponding
        VG
    threads : int
        number of threads to use during indexing
    verbose : bool
        print information about graph indexing
        
    Returns
    -------
    int
        status of VG indexing (0 = all ok; 1 = an error occurred)
    """

    if not isinstance(vg, str):
        raise ValueError("Invalid path to the genome variation graph. Exiting")

    if not os.path.exists(vg):
        errmsg = "unable to find {0}".format(vg)
        raise FileExistsError(errmsg)

    success: int

    # take chromosome name and add it the XG extension
    graph_name: str = vg.split('.')[-2]
    xg: str = ''.join([graph_name, ".xg"])
    gbwt: str = ''.join([graph_name, ".gbwt"])

    # perform indexing of the current genome variation graph
    if verbose:
        # print information about indexing 
        vg_index: str = 'vg index -t {0} -G {1} -v {2} -x {3} {4} -p'.format(threads, 
                                                                             gbwt, 
                                                                             vcf, 
                                                                             xg, 
                                                                             vg)
    else:
        vg_index = 'vg index -t {0} -G {1} -v {2} -x {3} {4}'.format(threads, gbwt, 
                                                                     vcf, xg, vg)

    code: int = subprocess.call(vg_index, shell=True)

    if code != 0:
        success = 1
    else:
        success = 0

    return success

# end of indexVG()


def get_chromlist(ref_genome: str) -> List[str]:
    """Scan the reference genome FASTA file to find the chromosomes for
    which there a sequence is available.
    
    The file must be in FASTA format and the chromosome names start with
    '>chr' (e.g. '>chrX', '>chr1', etc.)
        
    Parameters
    ----------
    ref_genome : str
        path to the reference genome FASTA file
        
    Returns
    -------
    list
        chomosomes for which a sequence is available in the given 
        reference genome FASTA file 
    """

    # redefine default SIGINT handler
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    # overwrite original SIGINT handler
    signal.signal(signal.SIGINT, original_sigint_handler)  
    chroms = list()

    print(
        "Reading the valid chromosome names from the given reference genome...\n"
        )

    try:
        with open(ref_genome, mode='r') as infile:
            for line in infile:
                line =  line.strip()
                if line[0] == ">":  # this line contains the chromosome name
                    if line[:4] == ">chr":
                        chroms.append(line[4:])  # remove the starting '>chr'
                    else:
                        chroms.append(line[1:])  # remove the starting '>

    except Exception as e:
        raise e 

    except KeyboardInterrupt:
        sigint_handler()

    finally:
        infile.close()  # close input stream

    return chroms

# end of get_chromslist()


def tbiexist(vcf: str) -> bool:
    """Check if the given VCF file has been indexed with tabix (chek the 
    existance of VCFname.vcf.gz.tbi file)

    Parameters
    ----------
    vcf : str
        path to the VCF file to be checked for a corresponding TBI index
        file 
        
    Returns
    -------
    bool
        TBI index search status
    
    """

    vcf_tbi: str = ''.join([vcf, '.tbi'])

    if os.path.isfile(vcf_tbi):
        return True

    return False

# end of tbiexist()

