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
from grafimo.utils import exception_handler, is_jaspar, is_meme, is_transfac, is_pfm
from grafimo.grafimo_errors import VGError, SubprocessError

import subprocess
import sys
import os


# GRAFIMO's version
__version__ = "1.1.6"


sys.excepthook = exception_handler  # custom exception handler


def buildvg(workflow: BuildVG, debug: bool) -> None:
    """Call the functions needed to constuct the genome variation graph 
    from a reference FASTA file and a phased VCF file.

    ...

    Parameters
    ----------
    workflow : BuildVG
        container of the argumentgs needed to build a genome variation
        graph

    Returns 
    -------
    None
    """

    if not isinstance(workflow, BuildVG):
        errmsg = f"Expected {type(BuildVG).__name__} object, got {type(workflow).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    printWelcomeMsg()
    # if verbose == True print a lot of info
    verbose = workflow.verbose
    print("\n\nBuilding the VG for chromosome:")
    for c in workflow.chroms:
        print(c, end=" ")
    print("\n")  # newline
    if verbose:
        print("Buildvg user parameters:")
        print("\t- Reference genome: ", workflow.reference_genome)
        print("\t- VCF file: ", workflow.vcf)
        print("\t- Reindex: ", workflow.reindex)
        print("\t- Chromosomes: ", workflow.chroms)
        print("\t- Chromosome prefix: ", workflow.chroms_prefix)
        print("\t- Name-map: ", workflow.namemap)
        print("\t- Cores: ", workflow.cores)
        print("\t- Output directory: ", workflow.outdir)
        print("\t- Debug:", debug)
        print("\t- Verbose: ", verbose)
        print("\t- Test mode: ", workflow.get_test())
    if verbose:
        print("\nBeginning VGs construction\n")
    # begin VGs construction
    construct_vg(workflow, debug)  # the VGs will be stored in the defined output directory

# end of buildvg()


def findmotif(workflow: Findmotif, debug: bool) -> None:
    """Call the functions needed to scan the genome variation graph for 
    the occurrence of a DNA motif, in genomic regions defined in the given 
    BED file.

    ...

    Parameters
    ----------
    workflow : Findmotif
        container of the arguments needed to scan the genome variation
        graph for the occurrences of the given DNA motif

    Returns
    -------
    None
    """
    
    if not isinstance(workflow, Findmotif):
        errmsg = f"Expected {type(Findmotif).__name__}, got {type(workflow).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    printWelcomeMsg()
    # if verbose == True print a lot of info
    verbose = workflow.verbose
    cores = workflow.cores  # cores to use
    if verbose:
        print("Findmotif user parameters:")
        if workflow.has_graphgenome():
            print("\t- Genome variation graph: ", workflow.graph_genome)
        elif workflow.has_graphgenome_dir():
            print("\t- Genome variation graphs directory: ", workflow.graph_genome_dir)
        else:  # both False
            errmsg = "Something went wrong during commandline args parse.\n"
            exception_handler(IOError, errmsg, debug)
        print("\t- BED file: ", workflow.bedfile)
        print("\t- Motif file: ", workflow.motif)
        print("\t- Chromosomes: ", workflow.chroms)
        print("\t- Chromome prefix: ", workflow.chroms_prefix)
        print("\t- Name-map: ", workflow.namemap)
        print("\t- Background: ", workflow.bgfile)
        print("\t- Pseudocount: ", workflow.pseudo)
        print("\t- Threshold: ", workflow.threshold)
        print("\t- Output directory: ", workflow.outdir)
        print("\t- Cores: ", cores)
        print("\t- recomb: ", workflow.recomb)
        print("\t- Top graphs: ", workflow.top_graphs)
        print("\t- no-qvalue: ", workflow.noqvalue)
        print("\t- no-reverse: ", workflow.noreverse)
        print("\t- Text only: ", workflow.text_only)
        print("\t- Threshold on q-values: ", workflow.qvalueT)
        print("\t- Verbose: ", verbose)
        print("\t- Test mode: ", workflow.get_test())
        print()  # newline
    # if genome graph is given, check that it is indexed (XG), otherwise index it
    if workflow.has_graphgenome():
        if workflow.graph_genome.split(".")[-1] == "vg":
            warnmsg = f"\nWARNING: {workflow.graph_genome} is not indexed. To scan a VG, it should be indexed."
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
                code = index_vg(workflow.graph_genome, vcf, cores, verbose, debug)  
                if code != 0:
                    errmsg = f"An error occurred during {workflow.graph_genome} indexing.\n"
                    exception_handler(VGError, errmsg, debug) 
            else:
                errmsg = f"To scan {workflow.graph_genome} are required the XG and GBWT indexes.\n"
                exception_handler(FileNotFoundError, errmsg, debug)
            vg_name = workflow.graph_genome.split(".vg")[-2]
            xg = ".".join([vg_name, "xg"])
            workflow.set_xg(xg)  # set indexed vg name (XG) 
    # motif PWM processing
    motifs = workflow.motif
    motif_set = MotifSet()
    for motif in motifs:
        if verbose:
            print(f"Processing motif(s) in: {motif}")
        m = get_motif_pwm(motif, workflow, cores, debug)
        if workflow.get_test():
            for i in range(len(m)):
                print(f"Score matrix of {motif}:")
                m[i].print("score_matrix")  
        motif_set.add_motif(m)
    # visit the graph to extract motif matches
    sequences_loc = scan_graph(motif_set.widths, workflow, debug)
    for motif in motif_set:
        # score sequences
        res = compute_results(motif, sequences_loc, debug, workflow)
        if workflow.text_only:
            print_results(res, debug)  # print results to stdout
        else:
            write_results(res, motif, motif_set.size, workflow, debug)
    # remove sequence directories
    if not workflow.get_test():
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

