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


from grafimo.utils import (
    sigint_handler, 
    exception_handler,
    ALL_CHROMS
)
from grafimo.workflow import BuildVG
from grafimo.grafimo_errors import VGError, SubprocessError, FileReadError

from typing import List, Dict

import subprocess
import signal
import time
import sys
import os


def get_reference_genome_from_ucsc(debug: bool) -> str:
    """Download the reference genome (hg38 assembly), from the UCSC
    database, in the current working directory and returns the path to 
    the corresponding FASTA file.

    ** WARNING **: This function has been written only for test purposes.

    ...
    
    Parameters
    ----------
    debug : bool
        Trace the complete error stack

    Returns
    -------
    str
        Path to the downloaded FASTA file (in .fa format)
    """

    # download genome
    address = "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    cmd = "wget -c {}".format(address)
    # the genome will be downloaded in the current directory
    code = subprocess.call(cmd, shell=True)  
    if code != 0:
        errmsg = "An error occurred while executing \"{}\". Exiting.\n"
        exception_handler(SubprocessError, errmsg.format(cmd), debug)
    # uncompress genome
    print("Uncompressing the genome...")
    genome_comp = "./hg38.fa.gz"
    if not os.path.exists(genome_comp):
        errmsg = f"Unable to find {genome_comp}.\n"
        exception_handler(FileNotFoundError, errmsg, debug)
    cmd = f"gunzip {genome_comp}"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        errmsg = f"An error occurred while executing \"{cmd}\". Exiting.\n"
        exception_handler(SubprocessError, errmsg, debug)
    # remove FASTA.GZ file if still present
    if os.path.exists(genome_comp):
        cmd = f"rm {genome_comp}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = f"An error occurred while executing \"{cmd}\". Exiting.\n"
            exception_handler(SubprocessError, errmsg, debug)
    # get the path to the genome file
    genome_uncomp = "./hg38.fa"  # should be in the current dir
    assert os.path.exists(genome_uncomp)
    genome = os.path.abspath(genome_uncomp)
    return genome

# end of get_reference_genome_from_ucsc()


def get_1000GProject_vcf(debug: bool) -> str:
    """Downloads a WGS VCF file from the 1000 Genome Project database
    (phase 3), containing SNVs and indels. The present file is used for 
    VG construction and graph indexing test purposes. 
    
    Since the variants present in this file are not phased, it cannot be 
    used to build the GBWT index and the corresponding haplotypes cannot
    be used. To use this features the VCF must be phased first.

    ** WARNING **: This function has been written only for test purposes.

    ...

    Parameters
    ----------
    debug : bool
        Trace the complete error stack
    
    Returns
    -------
    str
        Path to the downloaded VCF file (compressed)
    """

    # download the VCF
    address = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    address += "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"
    address += "ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz"
    cmd = f"wget -c {address}"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        errmsg = f"An error occurred while executing \"{cmd}\". Exiting.\n"
        exception_handler(SubprocessError, errmsg, debug)
    # vcf should be in the current dir
    vcf_file = "./ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz"
    assert os.path.exists(vcf_file)
    vcf = os.path.abspath(vcf_file)
    return vcf

# end of get1000GProject_vcf()


def construct_vg(buildvg_args: BuildVG, debug: bool) -> None:
    """ Create the genome graph using the given genome reference and 
    phased VCF file.
    
    The genome is not built as a single whole genome graph but a
    single graph is constructed for each chromosome.
    This approach avoids memory issues and allows the genome variation
    graph construction also on machines with less computational 
    resources.

    There is NO drawback using this approach with respect to building
    a whole genome graph and query it.
    
    Moreover, it allows parallel queries on the different chromosomes to
    be performed also on regular laptops (>= 16 GB of memory), which is 
    very difficult with a whole genome graph, that requires the user
    to set appropriately the number of cores to use. However, a whole
    genome graph can be still visited using a regular laptop using one 
    core.

    ...
    
    Parameters
    ----------
    buildvg_args : BuildVG
        Command line arguments required to construct the VG.
    debug : bool
        Trace the full error stack.

    Returns
    -------
    None 
    """
    
    if not isinstance(buildvg_args, BuildVG):
        errmsg = f"Expectd {type(BuildVG).__name__} object, got {type(buildvg_args).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    assert isinstance(debug, bool)
    # read the arguments to build the VGs
    reindex = buildvg_args.reindex
    chroms = buildvg_args.chroms
    chroms_prefix = buildvg_args.chroms_prefix
    namemap = buildvg_args.namemap
    threads = buildvg_args.cores
    outdir = buildvg_args.outdir
    verbose = buildvg_args.verbose
    test = buildvg_args.get_test()  # manually set in the code
    if test:
        reference = get_reference_genome_from_ucsc(debug)
        vcf = get_1000GProject_vcf(debug)
    else:
        reference = buildvg_args.reference_genome
        vcf = buildvg_args.vcf
    if verbose:
        print(f"using reference genome: {reference}")
        print(f"Using VCF file: {vcf}\n\n")
    if verbose:
        start_c = time.time()
        print(f"Reading chromosome names from {reference}...")
    # read chromosome names in reference FASTA
    chroms_available = get_chromlist(reference, debug)
    if verbose:
        end_c = time.time()
        print("done in %.2fs" % (end_c - start_c))
        print("Found chromosomes:\n", chroms_available, end="\n\n")
    if len(chroms) == 1 and chroms[0] == ALL_CHROMS:
        chroms = chroms_available
    else:
        # check user-defined chromosome names consistency with names in 
        # reference
        for c in chroms:
            if c not in chroms_available:
                errmsg = f"Chromosome \"{c}\" not found among names in {reference}.\n"
                exception_handler(ValueError, errmsg, debug)
    cwd = os.getcwd()  # recover current working directory
    # check if the VCF file has already been indexed with tabix
    if not tbiexist(vcf):
        msg = f"TBI file not found for {vcf.split('/')[-1]}. Indexing VCF file with tabix..."
        print(msg)
        cmd = f"tabix -p vcf {vcf}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:  # tabix didn't work
            errmsg = f"An error occurred while executing \"{cmd}\". Exiting.\n"
            exception_handler(SubprocessError, errmsg, debug)
    elif reindex:  # user asked to reindex VCF 
        msg = f"Reindexing {vcf.split('/')[-1]}...\n"
        print(msg)
        # remove old index
        cmd = f"rm {''.join([vcf, '.tbi'])}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = f"An error occurred while executing \"{cmd}\". Exiting.\n"
            exception_handler(SubprocessError, errmsg, debug)
        # reindex the VCF
        cmd = f"tabix -p vcf {vcf}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            # tabix didn't work
            errmsg = f"An error occurred while executing \"{cmd}\". Exiting.\n"
            exception_handler(SubprocessError, errmsg, debug)
    # enter the output directory
    os.chdir(outdir)
    if chroms_prefix: 
        assert not bool(namemap)
    if bool(namemap): 
        assert chroms_prefix != "chr"
    # build the VG for each chromosome or only for those told by user
    for chrname in chroms:
        if not bool(namemap): 
            chrom = "".join([chroms_prefix, chrname])
        elif bool(namemap): 
            try:
                chrom = namemap[chrname]
            except:
                errmsg = f"Missing out name map for chromosome \"{chrname}\".\n'"
                exception_handler(KeyError, errmsg, debug)
        vg = "".join([".", chrom, '.vg'])
        # build VG for current chromosome
        if verbose:
            start_build = time.time()
        code = build_vg(vg, reference, vcf, chrname, threads)
        if code != 0:
            errmsg = f"An error occurred during construction of {vg}.\n"
            exception_handler(VGError, errmsg, debug)
        if verbose:
            end_build = time.time()
            msg = f"Elapsed time to build {vg}:"   
            print(msg, "%.2fs" % (end_build - start_build), sep=" ")
        # index VG
        if verbose:
            start_index = time.time()
        msg = f"Indexing {chrom} VG and building the GBWT index..."
        print(msg)
        code = index_vg(vg, vcf, threads, verbose, debug)
        if code != 0:
            errmsg = f"An error occurred while indexing {vg}."
            exception_handler(VGError, errmsg, debug)
        if verbose:
            end_index = time.time()
            msg = f"Elapsed time to index {vg}"
            print(msg, "%.2fs" % (end_index - start_index), sep=" ")
        # The majority of applications work only with indexed graph,
        # so to save disk space is worth to delete the VGs and keep
        # only the XGs (is simple to go back using VG built-in functions)
        if verbose:
            print(f"Deleting {vg}")
        cmd = f"rm {vg}"
        subprocess.call(cmd, shell=True)
        if code != 0: 
            errmsg = f"An error occurred while executing \"{cmd}\". Exiting.\n"
            exception_handler(SubprocessError, errmsg, debug)
    # get VGs location
    graphs_loc = os.getcwd()
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

    ...
    
    Parameters
    ----------
    vg : str
        Path to the genome variation graph
    ref_genome : str
        Path to the reference genome FASTA file
    vcf : str
        Path to the phased VCF file (GZIP compressed)
    chrom_num : int
        Chromosomes for which GRAFIMO will build the VG.
    threads : int 
        Number of threads to use.
        
    Returns
    -------
    int  
        VG construction operation status (0 == ok; 1 == error)
    """

    build = f"vg construct -t {threads} -r {ref_genome} -v {vcf} -R {chrom_num} -C -a -p > {vg}"
    code = subprocess.call(build, shell=True)
    if code != 0:  # an error occurred during VG construction 
        success = 1
    else:  # construction OK
        success = 0
    return success

# end of build_vg()


def index_vg(
    vg: str, 
    vcf: str, 
    threads: int, 
    verbose: bool,
    debug: bool
) -> int:
    """Construct the XG and GBWT indexes for the given genome variation 
    graph. These indexes are required to query the genome when extracting
    motif occurrence candidates.

    The GBWT index allows to keep track of the haplotypes used to build
    the graph data structure and retrieve the samples genomes. 
    
    The indexing operation could take some time.

    ...
        
    Parameters
    ----------
    vg : str
        Path to the genome variation graph (VG format).
    vcf : str 
        Path to the phased VCF file.
    threads : int
        Number of threads to use.
    verbose : bool
        Print some information about graph indexing
    debug : bool
        Trace the full error stack.
        
    Returns
    -------
    int
        VG indexing operation status (0 == ok; 1 == error).
    """

    if not isinstance(vg, str):
        errmsg = f"Expected {str.__name__} instance, got {type(vg).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)       
    if not os.path.exists(vg):
        errmsg = f"Unable to find {vg}.\n"
        exception_handler(FileNotFoundError, errmsg, debug)

    # take chromosome name and add it the XG extension
    graph_name = vg.split('.')[-2]
    xg = "".join([graph_name, ".xg"])
    gbwt = "".join([graph_name, ".gbwt"])
    # perform indexing of the current genome variation graph
    if verbose:
        # print information about indexing 
        vg_index = f"vg index -t {threads} -G {gbwt} -v {vcf} -x {xg} {vg} -p"
    else:
        vg_index = f"vg index -t {threads} -G {gbwt} -v {vcf} -x {xg} {vg}"
    code = subprocess.call(vg_index, shell=True)
    if code != 0:
        success = 1
    else:
        success = 0
    return success

# end of index_vg()


def get_chromlist(ref_genome: str, debug: bool) -> List[str]:
    """Scan the reference genome FASTA file to find the chromosomes for
    which there a sequence is available.
    
    The file must be in FASTA format and the chromosome names start with
    '>chr' (e.g. '>chrX', '>chr1', etc.).

    ...
        
    Parameters
    ----------
    ref_genome : str
        Path to the reference genome FASTA file
    debug : bool
        Trace the full error stack.
        
    Returns
    -------
    List[str]
        Chromosomes found in the input FASTA file.
    """

    assert os.path.isfile(ref_genome)
    # redefine default SIGINT handler
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    # overwrite original SIGINT handler
    signal.signal(signal.SIGINT, original_sigint_handler)  
    chroms = list()
    try:
        with open(ref_genome, mode="r") as handle:
            while True:
                line = handle.readline()
                if not line: 
                    return  # empty file ?
                if line[0] == ">": 
                    break  # data start here
            while True:
                if line[0] != ">":
                    errmsg = "Sequence names in FASTA file should begin with \">\"\n."
                    exception_handler(FileReadError, errmsg, debug)
                else:
                    seqname = line.rstrip().split()[0][1:]  # skip ">"
                line = handle.readline()
                while True:
                    if not line: 
                        break  # empty sequence ?
                    if line[0] == ">": 
                        break  # sequence end
                    line = handle.readline()
                chroms.append(seqname)
                if not line: 
                    break  # reached EOF
    except KeyboardInterrupt:
        sigint_handler()
    except:
        errmsg = f"A problem was encountered reading {ref_genome}\n."
        exception_handler(FileReadError, errmsg, debug)
    finally:
        handle.close()
    return chroms

# end of get_chromslist()


def tbiexist(vcf: str) -> bool:
    """Check if already exists a TBI index for the input VCF file.

    ...

    Parameters
    ----------
    vcf : str
        Path to VCF file
        
    Returns
    -------
    bool
        Result
    """

    vcf_tbi = ".".join([vcf, "tbi"])
    if os.path.isfile(vcf_tbi):
        return True
    return False

# end of tbiexist()

