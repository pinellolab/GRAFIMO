"""From this point will be chosen the path to follow by GRAFIMO.

Here are called the functionalities to construct a genome variation 
graph from a reference genome FASTA file and a phased VCF file or to 
scan a precomputed graph for the occurrences of a given motif.

"""


from grafimo.workflow import BuildVG, Findmotif
from grafimo.constructVG import construct_vg, index_vg
from grafimo.motif_ops import get_motif_pwm
from grafimo.motif_set import MotifSet
from grafimo.extract_regions import scan_graph
from grafimo.res_writer import print_results, write_results
from grafimo.score_sequences import compute_results
from grafimo.utils import exception_handler
from grafimo.grafimo_errors import VGError, SubprocessError

import subprocess
import sys
import os


# GRAFIMO's version
__version__ = "1.1.5"


sys.excepthook = exception_handler  # custom exception handler


def buildvg(args_obj: BuildVG, debug: bool) -> None:
    """Call the functions needed to constuct the genome variation graph 
    from a reference FASTA file and a phased VCF file.

    ...

    Parameters
    ----------
    args_obj : BuildVG
        container of the argumentgs needed to build a genome variation
        graph

    Returns 
    -------
    None
    """

    if not isinstance(args_obj, BuildVG):
        errmsg = f"Expected {type(BuildVG).__name__} object, got {type(args_obj).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    printWelcomeMsg()
    # if verbose == True print a lot of info
    verbose = args_obj.verbose
    print("\n\nBuilding the VG for chromosome:")
    for c in args_obj.chroms:
        print(c, end=" ")
    print("\n")  # newline
    if verbose:
        print("Buildvg user parameters:")
        print("\t- Reference genome: ", args_obj.reference_genome)
        print("\t- VCF file: ", args_obj.vcf)
        print("\t- Reindex: ", args_obj.reindex)
        print("\t- Chromosomes: ", args_obj.chroms)
        print("\t- Chromosome prefix: ", args_obj.chroms_prefix)
        print("\t- Name-map: ", args_obj.namemap)
        print("\t- Cores: ", args_obj.cores)
        print("\t- Output directory: ", args_obj.outdir)
        print("\t- Debug:", debug)
        print("\t- Verbose: ", verbose)
        print("\t- Test mode: ", args_obj.get_test())
    if verbose:
        print("\nBeginning VGs construction\n")
    # begin VGs construction
    construct_vg(args_obj, debug)  # the VGs will be stored in the defined output directory

# end of buildvg()


def findmotif(args_obj: Findmotif, debug: bool) -> None:
    """Call the functions needed to scan the genome variation graph for 
    the occurrence of a DNA motif, in genomic regions defined in the given 
    BED file.

    ...

    Parameters
    ----------
    args_obj : Findmotif
        container of the arguments needed to scan the genome variation
        graph for the occurrences of the given DNA motif

    Returns
    -------
    None
    """
    
    if not isinstance(args_obj, Findmotif):
        errmsg = f"Expected {type(Findmotif).__name__}, got {type(args_obj).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    printWelcomeMsg()
    # if verbose == True print a lot of info
    verbose = args_obj.verbose
    cores = args_obj.cores  # cores to use
    if verbose:
        print("Findmotif user parameters:")
        if args_obj.has_graphgenome():
            print("\t- Genome variation graph: ", args_obj.graph_genome)
        elif args_obj.has_graphgenome_dir():
            print("\t- Genome variation graphs directory: ", args_obj.graph_genome_dir)
        else:  # both False
            errmsg = "Something went wrong during commandline args parse.\n"
            exception_handler(IOError, errmsg, debug)
        print("\t- BED file: ", args_obj.bedfile)
        print("\t- Motif file: ", args_obj.motif)
        print("\t- Chromosomes: ", args_obj.chroms)
        print("\t- Chromome prefix: ", args_obj.chroms_prefix)
        print("\t- Name-map: ", args_obj.namemap)
        print("\t- Background: ", args_obj.bgfile)
        print("\t- Pseudocount: ", args_obj.pseudo)
        print("\t- Threshold: ", args_obj.threshold)
        print("\t- Output directory: ", args_obj.outdir)
        print("\t- Cores: ", cores)
        print("\t- recomb: ", args_obj.recomb)
        print("\t- Top graphs: ", args_obj.top_graphs)
        print("\t- no-qvalue: ", args_obj.noqvalue)
        print("\t- no-reverse: ", args_obj.noreverse)
        print("\t- Text only: ", args_obj.text_only)
        print("\t- Threshold on q-values: ", args_obj.qvalueT)
        print("\t- Verbose: ", verbose)
        print("\t- Test mode: ", args_obj.get_test())
        print()  # newline
    # if genome graph is given, check that it is indexed (XG), otherwise index it
    if args_obj.has_graphgenome():
        if args_obj.graph_genome.split(".")[-1] == "vg":
            warnmsg = f"\nWARNING: {args_obj.graph_genome} is not indexed. To scan a VG, it should be indexed."
            # quick IO with user
            answer = "-1"  # ensure we enter while loop
            while answer.upper() != "Y" and answer.upper() != "N":
                print(warnmsg)
                answer = input("Do you want to index it now? (y/n)\n")
            if answer.upper() == "Y":
                vcf = input(
                    "Enter the path to the VCF file for haplotype tracking:\n"
                )
                if vcf.split(".")[-2] != "vcf" or vcf.split(".")[-1] != "gz":
                    errmsg = f"Unable to use {vcf} for haplotype tracking. The VCF must be compressed (e.g.\"myvcf.vcf.gz\".\n"
                    exception_handler(ValueError, errmsg, debug)
                if not os.path.isfile(vcf):
                    errmsg = f"Unable to locate {vcf}.\n"
                    exception_handler(FileNotFoundError, errmsg, debug)
                # VG indexing
                code = index_vg(args_obj.graph_genome, vcf, cores, verbose, debug)  
                if code != 0:
                    errmsg = f"An error occurred during {args_obj.graph_genome} indexing.\n"
                    exception_handler(VGError, errmsg, debug) 
            else:
                errmsg = f"To scan {args_obj.graph_genome} are required the XG and GBWT indexes.\n"
                exception_handler(FileNotFoundError, errmsg, debug)
            vg_name = args_obj.graph_genome.split(".vg")[-2]
            xg = ".".join([vg_name, "xg"])
            args_obj.set_xg(xg)  # set indexed vg name (XG) 
    # motif PWM processing
    motifs = args_obj.motif
    motif_set = MotifSet()
    for motif in motifs:
        if verbose:
            print(f"Processing motif(s) in: {motif}")
        m = get_motif_pwm(motif, args_obj, cores, debug)
        if args_obj.get_test():
            for i in range(len(m)):
                print(f"Score matrix of {motif}:")
                m[i].print("score_matrix")  
        motif_set.add_motif(m)
    # visit the graph to extract motif matches
    sequences_loc = scan_graph(motif_set.widths, args_obj, debug)
    for motif in motif_set:
        # score sequences
        res = compute_results(motif, sequences_loc, debug, args_obj)
        if args_obj.text_only:
            print_results(res, debug)  # print results to stdout
        else:
            write_results(res, motif, motif_set.size, args_obj, debug)
    # remove sequence directories
    if not args_obj.get_test():
        cmd = f"rm -rf {sequences_loc}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = f"An error occurred while executing {cmd}.\n"
            exception_handler(SubprocessError, errmsg, debug)

# end of findmotif()


def printWelcomeMsg() -> None:
    """Print GRAFIMO welcome message.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    for _ in range(75): print('*', end='')
    msg = "\n\tWelcome to GRAFIMO v{}"
    print()  # newline
    print(msg.format(__version__))
    print()  # newline
    for _ in range(75): print('*', end='')
    print()  # newline

# end of printWelcomeMsg()

