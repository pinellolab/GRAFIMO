"""Functions to retrieve genomic sequences from genome variation graphs.
Each query will extract the k-mers (k == motif width) from the VG,
in the given set of genomic regions.

Th genomic regions are given as set of coordinates in a UCSC BED file.

The k-mers are extracted on both forward and reverse complement starnds (by
default) or only from forward strand.

By default, are kept in the analysis only those k-mers appearing in the alternattive
genomes endoded in the scanned VG(s). Those k-mers follow the samples haplotypes.

If required, can also be extracted  all the possible recombinants k-mers, which 
can be obtained from the set of genetic variants encoded in the VG (--recomb) 
option. In this case the sample haplotypes are ignored.
"""


from grafimo.GRAFIMOException import SubprocessError, NotValidFFException, \
    FileReadError, VGError, FileFormatError
from grafimo.utils import die, ALL_CHROMS, printProgressBar, sigint_handler, exception_handler, NOMAP, isbed
from grafimo.score_sequences import ResultTmp
from grafimo.workflow import Findmotif
from grafimo.motif import Motif

from typing import List, Optional, Dict, Tuple, Set

import multiprocessing as mp

import subprocess
import warnings
import tempfile
import signal
import gzip
import time
import sys
import os


verbose = False  # global variable


def scan_graph(widths: Set[int], args_obj: Findmotif, debug: bool) -> str:
    if not isinstance(widths, set):
        errmsg = f"Expected {set.__name__}, got {type(widths).__name__}."
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(args_obj, Findmotif):
        errmsg = f"Expected {type(Findmotif).__name__}, got {type(args_obj).__name__}."
        exception_handler(TypeError, errmsg, debug)
    # recover arguments to begin VG scan
    if args_obj.has_graphgenome():  # single VG
        vg = args_obj.graph_genome
        if not is_vg_indexed(vg, debug):
            errmsg = "The genome variation graph should be indexed (XG format)."
            exception_handler(VGError, errmsg, debug)
    elif args_obj.has_graphgenome_dir():
        vg = args_obj.graph_genome_dir
    else:
        errmsg = "Unexpected genome variation graph found."
        exception_handler(VGError, errmsg, debug)
    bedfile = args_obj.bedfile
    chroms = args_obj.chroms
    chroms_prefix = args_obj.chroms_prefix
    namemap = args_obj.namemap
    cores = args_obj.cores
    global verbose  # change global verbose var value
    verbose = args_obj.verbose
    # begin sequence extraction
    try:
        print(f"\nExtracting regions defined in {bedfile}.\n")
        if verbose:
            start_bp = time.time()
        regions, region_num = get_regions_bed(bedfile, debug)
        if verbose:
            end_bp = time.time()
            print(
                "%s parsed in %.2fs. Found %d regions.\n" % (
                    bedfile, (end_bp - start_bp), region_num
                )
            )
        if args_obj.chroms_num == 1 and chroms[0] == ALL_CHROMS:
            chroms = [c.split("chr")[1] for c in regions.keys()]
        tmpwd = tempfile.mkdtemp(prefix="grafimo_") # create tmp directory
        cwd = os.getcwd()
        # enter the previously created temporary directory
        os.chdir(tmpwd)
        # create the queries list
        queries = []
        # overwrite the default SIGINT handler to exit gracefully
        # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool = mp.Pool(processes=cores)  # use #cores processes
        signal.signal(signal.SIGINT, original_sigint_handler)
        # create queries for each motif width
        for width in widths:
            # create directory with each motif width
            width_dir = f"width_{width}"
            os.mkdir(width_dir)
            if args_obj.has_graphgenome_dir():  # set of genome graphs
                for chrom in chroms:
                    if not bool(namemap):
                        chrname = "".join([chroms_prefix, chrom])
                    else:
                        try:
                            if chrom.startswith("chr"):
                                chrname = namemap[chrom.split("chr")[1]]
                            else:
                                chrname = namemap[chrom]
                        except:
                            errmsg = f"Missing name map for chromosome {chrom}.\n"
                            exception_handler(KeyError, errmsg, debug)
                    if chrom.startswith("chr"):
                        positions = regions[chrom]
                    else:
                        positions = regions["".join(["chr", chrom])]
                    for position in positions:
                        start = position[0]
                        stop = position[1]
                        if bool(namemap):
                            if chrom.startswith("chr"):
                                c = namemap[chrom.split("chr")[1]]
                            else:
                                c = chrom
                        elif chroms_prefix:
                            c = chrname.split(chroms_prefix)[1]
                        else:
                            c = chrname
                        region_index = "-".join([":".join([c, str(start)]), str(stop)])
                        region_name = region_index.replace(":", "_")
                        sequences = os.path.join(
                            ".", width_dir, ".".join(
                                [region_name, "tsv"]
                            )
                        )
                        xg = os.path.join(vg, ".".join([chrname, "xg"]))
                        # GBWT files must have the same prefix of corresponding XGs
                        gbwt = xg.replace("xg", "gbwt")
                        if not os.path.isfile(xg):
                            errmsg = f"Unable to locate {xg}. Are your VGs named with \"chr\"? Consider using --chroms-prefix-find or chroms-namemap-find.\n"
                            exception_handler(VGError, errmsg, debug)
                        if not os.path.isfile(gbwt):
                            errmsg = f"Unable to locate {gbwt}. Are your VGs named with \"chr\"? Consider using --chroms-prefix-find or chroms-namemap-find.\n"
                            exception_handler(VGError, errmsg, debug)
                        query = f"vg find -p {region_index} -x {xg} -H {gbwt} -K {width} -E > {sequences}"
                        queries.append(query)
            elif args_obj.has_graphgenome():  # single genome graph
                for chrom in chroms:
                    if not bool(namemap):
                        chrname = "".join([chroms_prefix, chrom])
                    else:
                        try:
                            chrname = namemap[chrom]
                        except:
                            errmsg = f"Missing name map for chromosome {chrom}.\n"
                            exception_handler(KeyError, errmsg, debug)
                    if chrom.startswith("chr"): 
                        if chrom not in regions.keys():
                            errmsg = f"{chrom} does not appear among the chromosomes available in {bedfile}.\n"
                            exception_handler(KeyError, errmsg, debug)
                        positions = regions[chrom]
                    else:
                        if ("".join(["chr", chrom])) not in regions.keys():
                            errmsg = f"{chrom} does not appear among the chromosomes available in {bedfile}.\n"
                            exception_handler(KeyError, errmsg, debug) 
                        positions = regions["".join(["chr", chrom])]
                    for position in positions:
                        start = position[0]
                        stop = position[1]
                        if chroms_prefix: 
                            c = chrname.split(chroms_prefix)[1]
                        else: 
                            c = chrname
                        region_index = "-".join([":".join([c, str(start)]), str(stop)])
                        region_name = region_index.replace(":", "_")
                        sequences = os.path.join(
                            ".", width_dir, ".".join(
                                [region_name, "tsv"]
                            )
                        )
                        xg = vg
                        # the GBWT must have the same prefix as XG
                        gbwt = xg.replace("xg", "gbwt") 
                        if not os.path.exists(xg):
                            errmsg = f"Unable to locate {xg}. Are your VGs named with \"chr\"? Consider using --chroms-prefix-find or chroms-namemap-find.\n"
                            exception_handler(VGError, errmsg, debug)
                        if not os.path.isfile(gbwt):
                            errmsg = f"Unable to locate {gbwt}. Are your VGs named with \"chr\"? Consider using --chroms-prefix-find or chroms-namemap-find.\n"
                            exception_handler(VGError, errmsg, debug)
                        query = f"vg find -p {region_index} -x {xg} -H {gbwt} -K {width} -E > {sequences}"
                        queries.append(query)
        get_kmers(queries, pool, debug, verbose)
    except:
        errmsg = f"An error occurred while scanning {vg}.\n"
        if args_obj.has_graphgenome_dir() or args_obj.has_graphgenome(): 
            exception_handler(VGError, errmsg, debug)
        else:
            errmsg = "Chromosome name mismatch. Check chromosome name consistency.\n"
            exception_handler(VGError, errmsg, debug)
    sequence_loc = os.getcwd()  
    os.chdir(cwd) 
    return sequence_loc

# end of scan_graph()



    
    








# def scan_graph(
#     motif: Motif,
#     args_obj: Findmotif,
#     debug: bool
# ) -> str:
#     """Obtain all the sequences of length K from the genome variation graph. 
#     K is the motif width.

#     The k-mers are extracted from the genomic regions defined in a UCSC BED file
#     or ENCODE narrowPeak file.

#     By default are extracted only those k-mers found on the alterantive genome
#     sequences encoded in the scanned genome variation graph(s). It is possible
#     to consider all the possible recombinant which can be obtained from the set 
#     of genetic variants encoded in the VG (--recomb option).

#     To perform k-mers extraction are followed the paths (haplotypes) encoded in 
#     VGs (defined as (V,E,P), where V are set of nodes, E the set of edges, and
#     P the set of paths or the haplotypes).

#     ...
    
#     Parameters
#     ----------
#     motif : Motif 
#         DNA motif
#     args_obj : Findmotif  
#         commandline arguments container
#     debug : bool
#         trace the full error stack
        
#     Returns
#     -------
#     str 
#         location of sequences files
#     """

#     if not isinstance(motif, Motif):
#         errmsg = "Expected Motif, got {}.\n"
#         exception_handler(TypeError, errmsg.format(type(motif).__name__), debug)
#     if not isinstance(args_obj, Findmotif):
#         errmsg = "Expected Findmotif, got {}.\n"
#         exception_handler(TypeError, errmsg.format(type(args_obj).__name__), debug)

#     if args_obj.has_graphgenome():  # single VG
#         vg = args_obj.graph_genome
#         if not is_vg_indexed(vg, debug):
#             errmsg = "The genome variation graph is not indexed, index it before proceeding.\n"
#             exception_handler(VGError, errmsg, debug)
#     elif args_obj.has_graphgenome_dir():
#         vg = args_obj.graph_genome_dir
#     else:
#         errmsg = "Unexpected genome variation graph given.\n"
#         exception_handler(VGError, errmsg, debug)

#     bedfile: str = args_obj.bedfile
#     chroms: List[str] =  args_obj.chroms
#     chroms_prefix: str = args_obj.chroms_prefix
#     namemap: dict = args_obj.namemap
#     cores: int = args_obj.cores
#     motif_width: int = motif.width
#     # modify global var value
#     global verbose
#     verbose = args_obj.verbose

#     # sequence extraction begin
#     try:
#         print("\nExtracting regions defined in {}.\n".format(bedfile))
#         if verbose: start_bp = time.time()
#         regions, region_num = get_regions_bed(bedfile, debug)
#         if verbose: 
#             end_bp = time.time()
#             print("%s parsed in %.2fs. Found %d regions.\n" % (bedfile, (end_bp - start_bp), region_num))
#         if args_obj.chroms_num == 1 and chroms[0] == ALL_CHROMS:
#             chroms = [c.split("chr")[1] for c in regions.keys()]
#         tmpwd: str = tempfile.mkdtemp(prefix='grafimo_')  # create a tmp dir
#         cwd: str = os.getcwd()  # get the current location 
#         os.chdir(tmpwd)  # enter the tmp dir 
#         # create a list of queries
#         queries: List[str] = list()  
#         # overwrite the default SIGINT handler to exit gracefully
#         # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
#         original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
#         pool: mp.Pool = mp.Pool(processes=cores)  # use no. cores processes
#         signal.signal(signal.SIGINT, original_sigint_handler) 
#         if args_obj.has_graphgenome_dir(): 
#             for chrom in chroms:
#                 if not bool(namemap):
#                     chrname = "".join([chroms_prefix, chrom])
#                 else:
#                     try:
#                         if chrom.startswith("chr"): 
#                             chrname = namemap[chrom.split("chr")[1]]
#                         else:
#                             chrname = namemap[chrom]
#                     except:
#                         errmsg = "Missing out name map for chromosome {}.\n"
#                         exception_handler(KeyError, errmsg.format(chrom), debug)
#                 if chrom.startswith("chr"): positions = regions[chrom]
#                 else: positions = regions["".join(["chr", chrom])]
#                 for pos in positions:
#                     start: int = pos[0]
#                     stop: int = pos[1]
#                     if bool(namemap):
#                         if chrom.startswith("chr"): c = namemap[chrom.split("chr")[1]]
#                         else: c = chrom
#                     elif chroms_prefix: c = chrname.split(chroms_prefix)[1]
#                     else: c = chrname
#                     region_index:str = "-".join(
#                         [":".join([c, str(start)]), str(stop)]
#                     )
#                     region_name: str = "-".join(
#                         ["_".join([chrname, str(start)]), str(stop)]
#                     )
#                     seqs: str = os.path.join(".", ".".join([region_name, "tsv"]))
#                     xg: str = os.path.join(vg, ".".join([chrname, "xg"]))
#                     # the GBWT must have the same prefix as XG
#                     gbwt: str = os.path.join(vg, ".".join([chrname, "gbwt"]))  
#                     if not os.path.isfile(xg):
#                         errmsg = "Unable to locate {}. Are your VGs named with \"chr\"? Consider using --chroms-prefix-find or chroms-namemap-find.\n"
#                         exception_handler(VGError, errmsg.format(xg), debug)
#                     if not os.path.isfile(gbwt):
#                         errmsg = "Unable to locate {}. Are your VGs named with \"chr\"? Consider using --chroms-prefix-find or chroms-namemap-find.\n"
#                         exception_handler(VGError, errmsg.format(gbwt), debug)
#                     query: str = "vg find -p {} -x {} -H {} -K {} -E > {}".format(
#                         region_index, xg, gbwt, motif_width, seqs
#                     )
#                     queries.append(query)
#             get_kmers(queries, pool, debug, verbose)

#         elif args_obj.has_graphgenome():
#             for chrom in chroms:
#                 if not bool(namemap):
#                     chrname = "".join([chroms_prefix, chrom])
#                 else:
#                     try:
#                         chrname = namemap[chrom]
#                     except:
#                         errmsg = "Missing out name map for chromosome {}.\n"
#                         exception_handler(KeyError, errmsg.format(chrom), debug)
#                 if chrom.startswith("chr"): 
#                     if chrom not in regions.keys():
#                         errmsg = "{} does not appear among the chromosomes available in {}.\n"
#                         exception_handler(KeyError, errmsg.format(chrom, bedfile), debug)
#                     positions = regions[chrom]
#                 else:
#                     if ("".join(["chr", chrom])) not in regions.keys():
#                         errmsg = "{} does not appear among the chromosomes available in {}.\n"
#                         exception_handler(KeyError, errmsg.format(chrom, bedfile), debug) 
#                     positions = regions["".join(["chr", chrom])]
#                 for pos in positions:
#                     start: int = pos[0]
#                     stop: int = pos[1]
#                     if chroms_prefix: c = chrname.split(chroms_prefix)[1]
#                     else: c = chrname
#                     region_index:str = "-".join(
#                         [":".join([c, str(start)]), str(stop)]
#                     )
#                     region_name: str = "-".join(
#                         ["_".join([chrname, str(start)]), str(stop)]
#                     )
#                     seqs: str = os.path.join(".", ".".join([region_name, "tsv"]))
#                     xg: str = vg
#                     xg_prefix: str = xg.split(".xg")[0]
#                     # the GBWT must have the same prefix as XG
#                     gbwt: str = ".".join([xg_prefix, "gbwt"]) 
#                     if not os.path.exists(xg):
#                         errmsg = "Unable to locate {}. Are your VGs named with \"chr\"? Consider using --chroms-prefix-find or chroms-namemap-find.\n"
#                         exception_handler(VGError, errmsg.format(xg), debug)
#                     if not os.path.isfile(gbwt):
#                         errmsg = "Unable to locate {}. Are your VGs named with \"chr\"? Consider using --chroms-prefix-find or chroms-namemap-find.\n"
#                         exception_handler(VGError, errmsg.format(gbwt), debug)
#                     query: str = "vg find -p {} -x {} -H {} -K {} -E > {}".format(
#                         region_index, xg, gbwt, motif_width, seqs
#                     )
#                     queries.append(query)
#             get_kmers(queries, pool, verbose)
#     except:
#         errmsg = "An error occurred while scanning {}.\n"
#         if args_obj.has_graphgenome_dir(): 
#             exception_handler(VGError, errmsg.format(xg), debug)
#         elif args_obj.has_graphgenome(): 
#             exception_handler(VGError, errmsg.format(vg), debug)
#         else:
#             errmsg = "Chromosome name mismatch. Check chromosome name consistency.\n"
#             exception_handler(VGError, errmsg, debug)
#     sequence_loc: str = os.getcwd()  
#     os.chdir(cwd) 

#     return sequence_loc

# # end of scan_graph()


























def get_kmers(
    queries: List[str], 
    pool: mp.Pool, 
    debug: bool,
    verbose: Optional[bool] = False,
) -> None:
    """Retrieve sequences from genome variation graph(s). The k-mers search is
    made in parallel creating #cores processes.

    ...

    Parameters
    ----------
    queries : list
        list of queries
    pool : multiprocessing.Pool
        pool ps
    debug : bool
        trace the full error stack
    verbose : bool, optional
        print additional information
    """

    if not isinstance(queries, list):
        errmsg = "Expected list, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(queries).__name__), debug)

    if verbose: start_re: float = time.time()
    try:
        res: mp.pool.MapResult = (pool.map_async(get_seqs, queries))
        if not verbose:
            it: int = 0
            while (True):
                if res.ready():
                    printProgressBar(
                        tot, tot, prefix='Progress:', suffix='Complete', length=50
                    )
                    break
                if it == 0: tot = res._number_left
                remaining = res._number_left
                printProgressBar(
                    (tot - remaining), tot, prefix='Progress:', suffix='Complete', length=50
                )
                time.sleep(1)
                it += 1
        ret: list = res.get(60 * 60 * 60)  # does not ignore signals
    except KeyboardInterrupt:
        pool.terminate()
        sigint_handler()
    else:
        pool.close()
        if verbose:
            end_re: float = time.time()
            print("Extracted sequences from all regions in %.2fs" % (end_re - start_re))

# end of get_kmers()


def get_seqs(query: str) -> None:
    """Retrieve k-mers within the genomic region. K is the motif width.
        
    Parameters
    ----------
    query : str
        region query
    """

    if verbose:
        region: str = query.split('-p')[1].split('-K')[0]
        print("Extracting sequences from region:", region)
    code: int = subprocess.call(query, shell=True)  
    if verbose:
        if code != 0:
            warnmsg: str = "A problem occurred while retrieving sequences from {}.\n".format(region)
            warnings.warn(warnmsg)

# end of get_seqs()


def is_vg_indexed(vg: str, debug: bool) -> bool:
    """Check if the genome variation graph has been indexed (XG format).

    ...

    Parameters
    ----------
    vg : str 
        path to genome variation graph
    debug : bool
        trace the full error stack
        
    Returns
    -------
    bool
        check result
    """

    if not isinstance(vg, str):
        errmsg = "Expected str, got{}.\n"
        exception_handler(TypeError, errmsg.format(type(vg).__name__), debug)
    if not os.path.isfile(vg):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(vg), debug)

    ff = vg.split(".")[-1]
    if ff == "xg": return True
    elif ff == "vg": return False
    else:  # unknown genome variation graph format
        errmsg = "Unknown genome variation graph format (VG or XG allowed).\n"
        exception_handler(VGError, errmsg, debug)

# end of is_vg_indexed()


def get_regions_bed(bedfile: str, debug: bool) -> Tuple[Dict, int]:
    """Read BED file and store genomic regions in a dictionary with the 
    chromosome numbers as keys. This allows to optimize VG cache loading.

    ...

    Parameters
    ----------
    bedfile : str 
        path to BED file
    debug : bool
        trace the full error stack
        
    Returns
    -------
    dict
        genomic regions grouped by chromosome
    int 
        number of genomic regions
    """

    if not isinstance(bedfile, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(bedfile).__name__), debug)
    if not os.path.isfile(bedfile):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(bedfile), debug)
    if not isbed(bedfile, debug):
        errmsg = "{} is not a UCSC BED file.\n"
        exception_handler(FileFormatError, errmsg.format(bedfile), debug)
    if os.stat(bedfile).st_size == 0:
        errmsg = "{} is empty.\n"
        exception_handler(FileReadError, errmsg.format(bedfile), debug)

    regions: Dict = dict()
    region_num: int = 0 
    gzipped = False 
    ff = bedfile.split(".")[-1]
    if ff == "gz": gzipped = True  # file is compressed
    try:
        if gzipped: ifstream = gzip.open(bedfile, mode="rt")
        else: ifstream = open(bedfile, mode="r")
        while True:
            line = ifstream.readline()
            if not line: break  # EOF or empty line?
            if line.startswith("chr"):  # data
                chrom, start, stop = line.strip().split()[:3]
                if chrom not in regions.keys():
                    regions.update({chrom:[(start, stop)]})
                else:
                    regions[chrom].append((start, stop))
                region_num += 1
    except:
        errmsg = "An error occurred while reading {}.\n"
        exception_handler(FileReadError, errmsg.format(bedfile), debug)
    finally:
        ifstream.close()
    
    return regions, region_num

# end of get_regions_bed() 

