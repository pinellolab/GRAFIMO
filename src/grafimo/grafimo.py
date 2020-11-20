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
import pandas as pd
import time
import sys
import os


# version of GRAFIMO
__version__ = '1.1.2'


def buildvg(args_obj: BuildVG) -> None:
    """Call the functions needed to constuct the genome variation graph 
    from a reference FASTA file and a phased VCF file.

    Parameters
    ----------
    args_obj : BuildVG
        container of the argumentgs needed to build a genome variation
        graph

    """

    if not isinstance(args_obj, BuildVG):
        raise ValueError("Unknown arguments object type")

    printWelcomeMsg()

    # if verbose == True will be printed a lot of unusefull stuff
    verbose = args_obj.get_verbose()

    print("\n\nBuilding the VG for chromosome:")
    for c in args_obj.get_chroms():
        print(c, end=" ")
    print("\n")  # newline

    if verbose:
        print("User parameters:")
        print("\t- Reference genome: ", args_obj.get_reference_genome())
        print("\t- VCF file: ", args_obj.get_vcf())
        print("\t- Reindex: ", args_obj.get_reindex())
        print("\t- Chromosomes: ", args_obj.get_chroms())
        print("\t- Cores: ", args_obj.get_cores())
        print("\t- Output directory: ", args_obj.get_outdir())
        print("\t- Verbose: ", verbose)
        print("\t- Test mode: ", args_obj.get_test())
    # end if

    if verbose:
        print("\nBeginning VGs construction\n")

    # begin VGs construction
    construct_vg(args_obj)  # the VGs will be stored in the defined output directory

# end of buildvg()


def findmotif(args_obj: Findmotif) -> None:
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
        raise ValueError("Unknown arguments object type")

    printWelcomeMsg()

    # if verbose == True will be printed a lot of additional information
    verbose: bool = args_obj.get_verbose()

    cores: int = args_obj.get_cores()  # cores to use

    if verbose:
        print("User parameters:")
        if args_obj.has_graph_genome():
            print("\t- Genome variation graph: ", args_obj.get_graph_genome())
        elif args_obj.has_graph_genome_dir():
            print("\t- Genome variation graphs directory: ", args_obj.get_graph_genome_dir())

        print("\t- BED file: ", args_obj.get_bedfile())
        print("\t- Motif file: ", args_obj.get_motif())
        print("\t- Chromosomes: ", args_obj.get_chroms())
        print("\t- Background: ", args_obj.get_bgfile())
        print("\t- Pseudocount: ", args_obj.get_pseudo())
        print("\t- Threshold: ", args_obj.get_threshold())
        print("\t- Output directory: ", args_obj.get_outdir())
        print("\t- Cores: ", cores)
        print("\t- recomb: ", args_obj.get_recomb())
        print("\t- Top graphs: ", args_obj.get_top_graphs())
        print("\t- no-qvalue: ", args_obj.get_no_qvalue())
        print("\t- no-reverse: ", args_obj.get_no_reverse())
        print("\t- Text only: ", args_obj.get_text_only())
        print("\t- Threshold on q-values: ", args_obj.get_qvalueT())
        print("\t- Verbose: ", verbose)
        print("\t- Test mode: ", args_obj.get_test())
        print()  # newline
    # end if

    errmsg: str

    # if we have a whole genome graph, check that it is indexed (XG)
    # index it, otherwise
    if args_obj.has_graph_genome():
        vg: str = args_obj.get_graph_genome()

        if vg.split('.')[-1] == 'vg':
            warnmsg: str = "\nWARNING: the given VG has not been indexed. "
            warnmsg += "To scan it, the VG must be indexed."
            print(warnmsg)
            answer = input("Do you want to index it now? (y/n) ")

            while answer.upper() != 'Y' and answer.upper() != 'N':
                print(warnmsg)
                answer: str = input("Do you want to index it now? (y/n) ")


            if answer.upper() == 'Y':
                vcf: str = input(
                    "Enter the path to the VCF file to use to include haplotypes:\n"
                    )
                
                if vcf.split(".")[-2] != "vcf" or vcf.split(".")[-1] != "gz":
                    errmsg = "\n\nERROR: the given file is not a zipped VCF file (*.vcf.gz)"
                    raise Exception(errmsg)
                
                if not os.path.isfile(vcf):
                    errmsg = "\n\nERROR: unable to find the given VCF file"
                    raise FileNotFoundError(errmsg)
                
                code: int = indexVG(vg, vcf, cores, verbose)

                if code != 0:
                    errmsg = ''.join(["\n\nERROR: unable to index ", vg, 
                                      ". Exiting"])
                    raise Exception(errmsg)
    
            else:
                errmsg = "\n\nTo proceed you need to provide a XG and a GBWT "
                errmsg += "index for your genome graph. GRAFIMO will exit"
                print(errmsg)
                sys.exit(1)
            
            # end if

            vg = vg.split('.vg')[-2]
            vg = ''.join([vg, '.xg'])
            args_obj.set_xg(vg)
        # end if
    # end if

    # process the motifs
    motifs: list = args_obj.get_motif()
    mtfSet: MotifSet = MotifSet()

    for mtf in motifs:
        if verbose:
            print("Processing motif(s) in: ", mtf, "\n")

        m = get_motif_pwm(mtf, args_obj, cores)

        if args_obj.get_test():
            for i in range(len(m)):
                print("Score matrix:")
                m[i].print("score_matrix")
        # end if

        mtfSet.addMotif(m)
    # end for

    motif_num: int = mtfSet.length()

    for mtf in mtfSet:
        # extract sequences
        sequence_loc: str = scan_graph(mtf, args_obj)

        # score sequences
        res: pd.DataFrame = compute_results(mtf, sequence_loc, args_obj)

        # write results
        if args_obj.get_text_only():
            print_results(res)
        else:
            write_results(res, mtf.getMotifID(), motif_num, args_obj)
        # end if
    # end for

# end of findmotif()


def printWelcomeMsg() -> None:
    """Prints the welcome message for GRAFIMO
  
    """
    for _ in range(50):
        print('*', end='')

    print()  # newline
    print("\n\tWELCOME TO GRAFIMO v", __version__, sep='')
    print()  # newline

    for _ in range(50):
        print('*', end='')
    print()  # newline

# end of printWelcomeMsg()

