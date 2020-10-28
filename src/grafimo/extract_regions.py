"""Functions to query a given genome variation graph or a set of graphs.
Each query will extract the k-mers (k == motif width) from the graph
in the given set of genomic regions.

Th genomic regions are given as set of coordinates in a UCSC BED file.

The k-mers are both from the forward and reverse strand of the genomic 
sequence.

By default, are extracted all the possible recombinant k-mers that can 
be obtained using the genomic variants defined in the used VCF file, in 
the queried regions. Then they are filtered to keep only those belonging 
to the haplotypes of the samples present in the VCF file used to build 
the genome variation graph.

The user can also decide to consider all possible recombinant using the
--recomb option, while calling GRAFIMO from the command line.
"""


from grafimo.GRAFIMOException import SubprocessError, NotValidFFException, \
    FileReadingException, VGException
from grafimo.utils import die, CHROMS_LIST, printProgressBar, sigint_handler
from grafimo.score_sequences import ResultTmp
from grafimo.workflow import Findmotif
from grafimo.motif import Motif
from typing import List, Optional, Dict, Tuple
import multiprocessing as mp
import subprocess
import signal
import warnings
import tempfile
import time
import sys
import os


verbose: bool = False  # global variable


def scan_graph(
    motif: Motif,
    args_obj: Findmotif
) -> str:
    """Obtain all the sequences of length K from the given genome
    variation graph, where K is the motif length.

    The sequences are obtained from the regions defined in the input
    BED file (UCSC BED file format).

    The sequences extracted correspond to all possible recombinant
    ones which can be obtained using the genomic variants given in the
    VCF file used to build the queried VG. 
    
    Then, they are filtered to keep only those beloning to haplotypes of 
    the samples from which the VCF file variants come from.

    The user can also decide to keep them all using the --recomb option. 
    
    Parameters
    ----------
    motif : Motif 
        DNA motif PWM to search on the VG
    args_obj : Findmotif  
        container of the arguments used during genome variation graph
        scanning
        
    Returns
    -------
    str 
        location of the files with the sequences extracted
    """

    errmsg: str
    vg: str
    chroms: List[str]

    # check the input arguments
    if not isinstance(motif, Motif):
        errmsg = "\n\nERROR: unknown motif object type"
        raise ValueError(errmsg)

    if not isinstance(args_obj, Findmotif):
        errmsg = "Unknown arguments object type. "
        errmsg += "Cannot scan the genome variation graph. Exiting"
        raise ValueError(errmsg)

    if args_obj.has_graph_genome():
        vg = args_obj.get_graph_genome()

        if not isGraph_genome_xg(vg):
            errmsg = "\n\nERROR: the genome variation graph is not in XG format"
            raise VGException(errmsg)
        # end if

    elif args_obj.has_graph_genome_dir():
        vg = args_obj.get_graph_genome_dir()

    else:
        raise VGException("\n\nERROR: the genome variation graph is missing")
    # end if

    bedfile: str = args_obj.get_bedfile()
    motif_width: int = motif.getWidth()
    cores: int = args_obj.get_cores()

    global verbose
    verbose = args_obj.get_verbose()

    print("\nExtracting regions defined in", bedfile, "\n")

    # read the regions where search the motif occurrences from the given 
    # BED file
    regions: Dict
    region_num: int
    regions, region_num = getBEDregions(bedfile)
    
    if(args_obj.get_chroms_num() == 1 and 
        args_obj.get_chroms()[0] == 'ALL_CHROMS'):
        chroms = list(regions.keys())
    else:
        chroms = [''.join(['chr', c]) for c in args_obj.get_chroms()]
    
    if verbose:
        print("\nFound", region_num, "regions in", bedfile)


    # create a tmp working directory
    tmpwd: str = tempfile.mkdtemp(prefix='grafimo_')

    # get the new location of graphs wrt the tmp dir
    cwd: str = os.getcwd()

    # enter the tmp dir where store the extracted sequences
    os.chdir(tmpwd)

    # list of queries
    queries: List[str] = list()  

    # redefine default SIGINT handler
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool: mp.Pool = mp.Pool(processes=cores)  # use no. cores processes
     # overwrite the default SIGINT handler to exit gracefully
    # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
    signal.signal(signal.SIGINT, original_sigint_handler) 

    positions: List[Tuple[int, int]]

    if args_obj.has_graph_genome_dir():

        # vg -> directory containing a set of VGs
        if vg[-1] == "/":
            pass
        else:
            vg = ''.join([vg, "/"])
        # end if

        for chrom in chroms:
            positions = regions[chrom]

            for pos in positions:
                start: int = pos[0]
                stop: int = pos[1]

                # the chromosome is among the ones to query
                region_index:str = ''.join([chrom.split('chr')[-1], ':', 
                                            str(start), '-', str(stop)])
                region_name: str = ''.join([chrom, '_', str(start), '-', 
                                            str(stop)])
                seqs: str = os.path.join('.', ''.join([region_name, '.tsv']))

                xg: str = ''.join([vg, chrom, '.xg'])
                # the GBWT must have the same prefix as XG
                gbwt: str = ''.join([vg, chrom, '.gbwt'])  

                if not os.path.exists(xg):
                    errmsg = ''.join(["\n\nERROR: unable to use ", xg, 
                                      ". Exiting"])
                    raise FileNotFoundError(errmsg)

                if not os.path.isfile(gbwt):
                    errmsg = "ERROR: unable to find GBWT file for"
                    errmsg = ' '.join([errmsg, xg])
                    raise FileNotFoundError(errmsg)

                query: str = 'vg find -p {0} -x {1} -H {2} -K {3} -E > {4}'.format(region_index, 
                                                                                   xg, 
                                                                                   gbwt, 
                                                                                   motif_width, 
                                                                                   seqs)
                queries.append(query)

        # extract from the graph the binding site candidates to score
        get_kmers(queries, pool, verbose)

    elif args_obj.has_graph_genome():

        for chrom in chroms:
            positions = regions[chrom]

            for pos in positions:
                start: int = pos[0]
                stop: int = pos[1]

                # the chromosome is among the ones to query
                region_index: str = ''.join([chrom.split('chr')[-1], ':', 
                                             str(start), '-', str(stop)])
                region_name: str = ''.join([chrom, '_', str(start), '-', 
                                            str(stop)])
                seqs: str = os.path.join('.', ''.join([region_name, '.tsv']))

                xg: str = vg
                xg_prefix: str = xg.split(".xg")[-2]
                # the GBWT must have the same prefix as XG
                gbwt: str = ''.join([xg_prefix, '.gbwt']) 

                if not os.path.exists(xg):
                    errmsg = ''.join(["\n\nERROR: unable to use ", xg, ". Exiting"])
                    raise FileNotFoundError(errmsg)

                if not os.path.isfile(gbwt):
                    errmsg = "ERROR: unable to find GBWT file for"
                    errmsg = ' '.join([errmsg, xg])
                    raise FileNotFoundError(errmsg)

                query = 'vg find -p {0} -x {1} -H {2} -K {3} -E > {4}'.format(region_index, 
                                                                              xg, 
                                                                              gbwt, 
                                                                              motif_width, 
                                                                              seqs)
                queries.append(query)

        # extract from the graph the binding site candidates to score
        get_kmers(queries, pool, verbose)

    else:
        raise Exception("\n\nERROR: do not know how to proceed. Exiting")
    # end if

    # the extracted sequences are store in the cwd
    sequence_loc: str = os.getcwd()  
    os.chdir(cwd) 

    return sequence_loc

# end of scan_graph()


def get_kmers(
    queries: List[str], 
    pool: mp.Pool,  
    verbose: Optional[bool] = False
) -> None:
    """Extract the genomic sequences (both from reverse and forward 
    strands)in the queried regions from the VG.

    The sequence extraction is perfromed in parallel working on a 
    user defined number of cores (by default all the cores available).

    Parameters
    ----------
    queries : list
        set of queries to perform on the graph to extract the motif
        occurrence candidates
    pool : multiprocessing.Pool
        pool of parallel processes to run  
    verbose : bool, optional
        flag used to define if additional information has to printed

    """

    if not isinstance(queries, list):
        raise Exception

    if verbose:
        start_re: float = time.time()

    # extract regions
    try:
        # query the VGs
        res: mp.pool.MapResult = (pool.map_async(get_seqs, queries))

        if not verbose:
            it: int = 0
            while (True):
                if res.ready():
                    # when finished call for the last time 
                    # printProgressBar()
                    printProgressBar(tot, tot, prefix='Progress:',
                                    suffix='Complete', length=50)
                    break
                # end if
                if it == 0:
                    tot = res._number_left

                remaining = res._number_left
                printProgressBar((tot - remaining), tot, prefix='Progress:',
                                  suffix='Complete', length=50)
                time.sleep(2)
                it += 1
            # end while
        # end if

        ret: list = res.get(60 * 60 * 60)  # does not ignore signals

    except KeyboardInterrupt:
        pool.terminate()
        sigint_handler()

    else:
        pool.close()

        if verbose:
            end_re: float = time.time()
            print(
                "Extracted sequences from all regions in %.2fs" % (end_re - start_re)
                )
        # end if
    # end try


def get_xg_loc(toXGpath: str) -> str:
    """Get the path to the whole genome variation graphXG index

    Parameters
    ----------
    toXGpath : str
        path to the genome variation graph XG index
        
    Returns
    -------
    str 
        path to the directory containing the graph XG index
    """

    for i in range(1, len(toXGpath)):
        if toXGpath[-i] == '/':
            bp = -i
            break
    # end for

    toXGpath_dir: str = toXGpath[:bp]

    return toXGpath_dir
# end of get_xg_loc()


def get_seqs(query: str) -> None:
    """Retrieve the k-mers withing the current genomic region, where k
    is the motif length.

    The sequences are retrieved by calling the vg find built-in method
    of VG.
        
    Parameters
    ----------
    query : str
        query to apply on the genome variation graph        

    """

    if verbose:
        region: str = query.split('-p')[1].split('-K')[0]
        print("Extracting sequences from region:", region)

    code: int = subprocess.call(query, shell=True)  # perform query

    if verbose:
        if code != 0:
            warnmsg: str = ''.join(
                ["A problem occurred during sequences extraction in ", region]
                )
            warnings.warn(warnmsg)
        # end if
    # end if

# end of get_seqs()


def isGraph_genome_xg(vg: str) -> bool:
    """Check if the given genome variation graph is in XG format.

    To have a genome variation graph in XG format it must have been 
    indexed.

    Parameters
    ----------
    vg : str 
        path to a genome variation graph
        
    Returns
    -------
    bool
        is the given genome variation graph in XG format
    """

    errmsg: str

    if not isinstance(vg, str):
        errmsg = "\n\nERROR: Invalid path to the genome graph. Cannot proceed"
        raise VGException(errmsg)

    if vg.split('.')[-1] == 'xg':
        return True
    elif vg.split('.') == 'vg':
        return False
    else:
        errmsg = "\n\nERROR: do not know what to do with the given genome graph." 
        errmsg += " Only XG or VG format allowed"
        raise VGException(errmsg)

# end of isGraph_genome_xg()


def getBEDregions(bedfile: str) -> Tuple[Dict, int]:
    """Read the BED file containing the genomic regions to scan for the 
    occurrences of the given motif.

    The regions are stored in Dict with the chromosome numbers as keys.
    This allows less overhead while loading the same genome variation
    graph on the cache when executing sequences extraction. 

    Parameters
    ----------
    bedfile : str 
        path to the BED file
        
    Returns
    -------
    dict
        regions defined in the BED file grouped by chromosome
    int 
        number of regions contained in the BED file
    """

    if bedfile.split('.')[-1] != 'bed':  # not a BED file
        raise NotValidFFException("The given file is not a BED file")

    regions:Dict
    region_num: int
    chrom: str
    start: int
    stop: int

    regions = dict()
    region_num = 0 

    try:
        with open(bedfile, mode='r') as inbed:  
            for line in inbed:
                bedline = line.split()
                
                if len(bedline) < 3 or len(bedline) > 12:
                    errmsg: str = '\n'.join(
                        ["The given BED file is not in UCSC BED file format.",
                         "Please refer to https://genome.ucsc.edu/FAQ/FAQformat.html#format1\n"]
                         )
                    raise NotValidFFException(errmsg)
                
                chrom, start, stop = bedline[0:3]
                if chrom not in regions.keys():
                    regions.update({chrom: [(start, stop)]})
                else:
                    regions[chrom].append((start, stop))

                region_num += 1
                # end if
            # end for
        # end with

    except:  # not able to read the BED file
        msg: str = ' '.join(["\n\nError: unable to read", bedfile])
        raise FileReadingException(msg)

    else:
        return regions, region_num

    finally:
        inbed.close()  # close the file stream
    # end try
# end of getBEDregions()

