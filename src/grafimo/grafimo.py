"""From this point will be chosen the path to follow by GRAFIMO.

Here are called the functionalities to construct a genome variation 
graph from a reference genome FASTA file and a phased VCF file or to 
scan a precomputed graph for the occurrences of a given motif.

"""


from grafimo.workflow import BuildVG, Findmotif
from grafimo.constructVG import construct_vg, indexVG
from grafimo.motif_ops import get_motif_pwm
from grafimo.motif_set import MotifSet
from grafimo.extract_regions import scan_graph
from grafimo.res_writer import print_results, write_results
from grafimo.score_sequences import compute_results
from grafimo.utils import exception_handler

import pandas as pd

import time
import sys
import os


# version of GRAFIMO
__version__ = '1.1.4'


sys.excepthook = exception_handler


def buildvg(args_obj: BuildVG, debug: bool) -> None:
    """Call the functions needed to constuct the genome variation graph 
    from a reference FASTA file and a phased VCF file.

    Parameters
    ----------
    args_obj : BuildVG
        container of the argumentgs needed to build a genome variation
        graph

    """

    if not isinstance(args_obj, BuildVG):
        errmsg = "Expected BuildVG object, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(args_obj).__name__), debug)
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
    # end if

    if verbose:
        print("\nBeginning VGs construction\n")

    # begin VGs construction
    construct_vg(args_obj, debug)  # the VGs will be stored in the defined output directory

# end of buildvg()


def findmotif(args_obj: Findmotif, debug: bool) -> None:
    """Call the functions needed to scan the genome variation graph for 
    the occurrence of a DNA motif, in genomic regions defined in the given 
    BED file.

    Parameters
    ----------
    args_obj : Findmotif
        container of the arguments needed to scan the genome variation
        graph for the occurrences of the given DNA motif

    """
    
    if not isinstance(args_obj, Findmotif):
        errmsg = "Expected Findmotif, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(args_obj).__name__), debug)
    printWelcomeMsg()
    # if verbose == True print a lot of info
    verbose: bool = args_obj.verbose
    cores: int = args_obj.cores  # cores to use

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
    # end if

    errmsg: str
    # if genome graph is given, check that it is indexed (XG), otherwise index it
    if args_obj.has_graphgenome():
        if args_obj.graph_genome.split('.')[-1] == 'vg':
            warnmsg: str = "\nWARNING: {} is not indexed. To scan a VG, it should be indexed."
            # quick IO with user
            answer = -1  # ensure we enter while loop
            while answer.upper() != 'Y' and answer.upper() != 'N':
                print(warnmsg.format(args_obj.graph_genome))
                answer: str = input("Do you want to index it now? (y/n)\n")
            if answer.upper() == 'Y':
                vcf: str = input(
                    "Enter the path to the VCF file for haplotype tracking:\n"
                )
                if vcf.split(".")[-2] != "vcf" or vcf.split(".")[-1] != "gz":
                    errmsg = "Unable to use {} for haplotype tracking. The VCF must be compressed (e.g.\"myvcf.vcf.gz\".\n"
                    exception_handler(ValueError, errmsg.format(vcf), debug)
                if not os.path.isfile(vcf):
                    errmsg = "Unable to locate {}.\n"
                    exception_handler(FileNotFoundError, errmsg.format(vcf), debug)
                code: int = indexVG(args_obj.graph_genome, vcf, cores, verbose)  # VG indexing
                if code != 0:
                    errmsg = "An error occurred during {} indexing.\n"
                    exception_handler(VGException, errmsg.format(args_obj.graph_genome), debug) 
            else:
                errmsg = "To scan {} are required the XG and GBWT indexes.\n"
                exception_handler(FileNotFoundError, errmsg.format(args_obj.graph_genome), debug)
            vg_name = args_obj.graph_genome.split('.vg')[-2]
            xg = ".".join([vg_name, "xg"])
            args_obj.set_xg(xg)  # set indexed vg name (XG) 
        # end if
    # end if

    # motif PWM processing
    motifs: list = args_obj.motif
    mtfSet: MotifSet = MotifSet()
    for mtf in motifs:
        if verbose:
            print("Processing motif(s) in: {}".format(mtf))
        m = get_motif_pwm(mtf, args_obj, cores, debug)
        if args_obj.get_test():
            for i in range(len(m)):
                print("Score matrix of {}:".format(mtf))
                m[i].print("score_matrix")  
        mtfSet.addMotif(m)
    for mtf in mtfSet:
        # extract sequences
        sequence_loc: str = scan_graph(mtf, args_obj, debug)
        # score sequences
        res: pd.DataFrame = compute_results(mtf, sequence_loc, debug, args_obj)
        # write results
        if args_obj.text_only: print_results(res)  # print to stdout
        else: write_results(res, mtf, mtfSet.size, args_obj, debug)

# end of findmotif()


def printWelcomeMsg() -> None:
    """Print GRAFIMO welcome message"""
    for _ in range(75): print('*', end='')
    msg = "\n\tWelcome to GRAFIMO v{}"
    print()  # newline
    print(msg.format(__version__))
    print()  # newline
    for _ in range(75): print('*', end='')
    print()  # newline

# end of printWelcomeMsg()

