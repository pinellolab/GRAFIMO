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


from grafimo.utils import die, sigint_handler, ALL_CHROMS, exception_handler
from grafimo.workflow import BuildVG
from grafimo.GRAFIMOException import VGError, ValueException, \
    SubprocessError, FileReadError
from typing import List, Dict
import subprocess
import signal
import time
import sys
import os


def get_reference_genome_from_ucsc(debug) -> str:
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
    address = "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    cmd = "wget -c {}".format(address)
    # the genome will be downloaded in the current directory
    code = subprocess.call(cmd, shell=True)  
    if code != 0:
        errmsg = "An error occurred while executing \"{}\". Exiting.\n"
        exception_handler(SubprocessError, errmsg.format(cmd), debug)
    # uncompress genome
    print("Uncompressing the genome...")
    genome_comp: str = './hg38.fa.gz'
    if not os.path.exists(genome_comp):
        errmsg = "Unable to find {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(genome_comp), debug)
    cmd = 'gunzip {0}'.format(genome_comp)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        errmsg = "An error occurred while executing \"{}\". Exiting.\n"
        exception_handler(SubprocessError, errmsg.format(cmd), debug)
    # remove FASTA.GZ file if still present
    if os.path.exists(genome_comp):
        cmd = 'rm {0}'.format(genome_comp)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = "An error occurred while executing \"{}\". Exiting.\n"
            exception_handler(SubprocessError, errmsg.format(cmd), debug)
    # get the path to the genome file
    genome_uncomp: str = "./hg38.fa"  # should be in the current dir
    assert os.path.exists(genome_uncomp)
    genome: str = os.path.abspath(genome_uncomp)

    return genome

# end of get_reference_genome_from_ucsc()


def get_1000GProject_vcf(debug) -> str:
    """Downloads a WGS VCF file from the 1000 Genome Project database
    (phase 3), containing SNVs and indels. The present file is used for 
    VG construction and graph indexing test purposes. 
    
    Since the variants present in this file are not phased, it cannot be 
    used to build the GBWT index and the corresponding haplotypes cannot
    be used. To use this features the VCF must be phased first.

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
        errmsg = "An error occurred while executing \"{}\". Exiting.\n"
        exception_handler(SubprocessError, errmsg.format(cmd), debug)
    # vcf should be in the current dir
    vcf_file: str = './ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz'
    assert os.path.exists(vcf_file)
    vcf: str = os.path.abspath(vcf_file)

    return vcf

# end of get1000GProject_vcf()


def construct_vg(buildvg_args: BuildVG, debug: bool) -> None:
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
        errmsg = "Expectd BuildVG object, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(buildvg_args).__name__), debug)
    
    # read the arguments to build the VGs
    reindex: bool = buildvg_args.reindex
    chroms: List[str] = buildvg_args.chroms
    chroms_prefix: str = buildvg_args.chroms_prefix
    namemap: Dict = buildvg_args.namemap
    threads: int = buildvg_args.cores
    outdir: str = buildvg_args.outdir
    verbose: bool = buildvg_args.verbose
    test: bool = buildvg_args.get_test()  # manually set in the code
    msg: str
    reference: str
    vcf: str

    if test:
        reference = get_reference_genome_from_ucsc(debug)
        vcf = get_1000GProject_vcf(debug)
    else:
        reference = buildvg_args.reference_genome
        vcf = buildvg_args.vcf

    if verbose:
        print("using reference genome: ", reference)
        print("Using VCF file: ", vcf, "\n\n")

    if verbose:
        start_c: float = time.time()
        print("Reading chromosome names from {}...".format(reference))
    # read chromosome names in reference FASTA
    chroms_available: List[str] = get_chromlist(reference, debug)
    if verbose:
        end_c: int = time.time()
        print("done in %.2fs" % (end_c - start_c))
        print("Found chromosomes:\n", chroms_available, end="\n\n")
    if len(chroms) == 1 and chroms[0] == ALL_CHROMS:
        chroms: List[str] = chroms_available
    else:
        # check user-defined chromosome names consistency with names in 
        # reference
        for c in chroms:
            if c not in chroms_available:
                errmsg = "Chromosome \"{}\" not found among names in {}.\n"
                exception_handler(ValueError, errmsg.format(c, reference), debug)

    cwd: str = os.getcwd()
    cmd: str
    code: int
    # check if the VCF file has already been indexed with tabix
    if not tbiexist(vcf):
        msg = "TBI file not found for {}. Indexing VCF file with tabix..."
        print(msg.format(vcf.split('/')[-1]))
        cmd = 'tabix -p vcf {0}'.format(vcf)
        code = subprocess.call(cmd, shell=True)
        if code != 0:  # tabix didn't work
            errmsg = "An error occurred while executing \"{}\". Exiting.\n"
            exception_handler(SubprocessError, errmsg.format(cmd), debug)
    elif reindex:  # user asked to reindex VCF 
        msg = "Reindexing {}...\n"
        print(msg.format(vcf.split('/')[-1]))
        # remove old index
        cmd = "rm {0}".format(''.join([vcf, ".tbi"]))
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = "An error occurred while executing \"{}\". Exiting.\n"
            exception_handler(SubprocessError, errmsg.format(cmd), debug)
        # reindex the VCF
        cmd = "tabix -p vcf {0}".format(vcf)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            # tabix didn't work
            errmsg = "An error occurred while executing \"{}\". Exiting.\n"
            exception_handler(SubprocessError, errmsg.format(cmd), debug)
    # end if

    # enter the output directory
    os.chdir(outdir)
    if chroms_prefix: assert not bool(namemap)
    if bool(namemap): assert chroms_prefix != "chr"
    # build the VG for each chromosome or only for those told by user
    for chrname in chroms:
        if not bool(namemap): 
            chrom: str = "".join([chroms_prefix, chrname])
        elif bool(namemap): 
            try:
                chrom: str = namemap[chrname]
            except:
                errmsg = "Missing out name map for chromosome \"{}\".\n'"
                exception_handler(KeyError, errmsg.format(chrname), debug)
        vg: str = ''.join([".", chrom, '.vg'])
        # build VG for current chromosome
        if verbose:
            start_build: float = time.time()
        code = build_vg(vg, reference, vcf, chrname, threads)
        if code != 0:
            errmsg = "An error occurred during construction of {}.\n"
            exception_handler(VGError, errmsg.format(vg), debug)
        if verbose:
            end_build: float = time.time()
            msg = "Elapsed time to build {}:"   
            print(msg.format(vg), "%.2fs" % (end_build - start_build), sep=" ")
        # index VG
        if verbose:
            start_index: float = time.time()
        msg = "Indexing {} VG and building the GBWT index..."
        print(msg.format(chrom))
        code = indexVG(vg, vcf, threads, verbose, debug)
        if code != 0:
            errmsg = "An error occurred while indexing {}."
            exception_handler(VGError, errmsg.format(vg), debug)
        if verbose:
            end_index: float = time.time()
            msg = "Elapsed time to index {}"
            print(msg.format(vg), "%.2fs" % (end_index - start_index), sep=" ")
        # end if

        # The majority of applications work only with indexed graph,
        # so to save disk space is worth to delete the VGs and keep
        # only the XGs (is simple to go back using VG built-in functions)
        if verbose:
            print("Deleting {0}".format(vg))
        cmd = 'rm {0}'.format(vg)
        subprocess.call(cmd, shell=True)
        if code != 0: 
            errmsg = "An error occurred while executing \"{}\". Exiting.\n"
            exception_handler(SubprocessError, errmsg.format(cmd), debug)
    # end for

    # get VGs location
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

    build = 'vg construct -t {4} -r {1} -v {2} -R {0} -C -a -p > {3}'.format(
        chrom_num, ref_genome, vcf,vg, threads
    )
    code: int = subprocess.call(build, shell=True)
    if code != 0:
        # an error occurred during VG construction 
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
    verbose: bool,
    debug: bool
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
        errmsg = "Expected str instance, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(vg).__name__), debug)       
    if not os.path.exists(vg):
        errmsg = "Unable to find {}.\n"
        exception_handler(FileNotFoundError, errmsg, debug)
    success: int
    # take chromosome name and add it the XG extension
    graph_name: str = vg.split('.')[-2]
    xg: str = ''.join([graph_name, ".xg"])
    gbwt: str = ''.join([graph_name, ".gbwt"])
    # perform indexing of the current genome variation graph
    if verbose:
        # print information about indexing 
        vg_index: str = 'vg index -t {0} -G {1} -v {2} -x {3} {4} -p'.format(
            threads, gbwt, vcf, xg, vg
        )
    else:
        vg_index = 'vg index -t {0} -G {1} -v {2} -x {3} {4}'.format(
            threads, gbwt, vcf, xg, vg
        )
    code: int = subprocess.call(vg_index, shell=True)
    if code != 0:
        success = 1
    else:
        success = 0

    return success

# end of indexVG()


def get_chromlist(ref_genome: str, debug: bool) -> List[str]:
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

    assert os.path.isfile(ref_genome)
    # redefine default SIGINT handler
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    # overwrite original SIGINT handler
    signal.signal(signal.SIGINT, original_sigint_handler)  
    chroms = list()
    
    try:
        with open(ref_genome, mode='r') as ifstream:
            while True:
                line = ifstream.readline()
                if not line: return  # empty file ?
                if line[0] == ">": break  # data start here
            while True:
                if line[0] != ">":
                    errmsg = "Sequence names in FASTA file should begin with \">\"\n."
                    exception_handler(FileReadError, errmsg, debug)
                else:
                    seqname = line.rstrip().split()[0][1:]  # skip ">"
                line = ifstream.readline()
                while True:
                    if not line: break  # empty sequence ?
                    if line[0] == ">": break  # sequence end
                    line = ifstream.readline()
                chroms.append(seqname)
                if not line: break  # reached EOF
    except KeyboardInterrupt:
        sigint_handler()
    except:
        errmsg = "A problem was encountered reading {}\n."
        exception_handler(FileReadError, errmsg.format(ref_genome), debug)
    finally:
        ifstream.close()

    return chroms

# end of get_chromslist()


def tbiexist(vcf: str) -> bool:
    """Check if already exists an index for the VCF file.

    Parameters
    ----------
    vcf : str
        path to VCF file
        
    Returns
    -------
    bool
        Check result
    
    """

    vcf_tbi: str = ".".join([vcf, "tbi"])
    if os.path.isfile(vcf_tbi):
        return True
    return False

# end of tbiexist()

