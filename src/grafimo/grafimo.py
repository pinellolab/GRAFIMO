"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it


Here are defined the steps of the two available pipelines:

    - with vg creation
    - only scoring (indexing of the vg if required [vg => xg])

The final results can be visualized directly on the terminal if
text_only == True (only in TSV format).
By default the results are written to the repository given by the user.

"""


from grafimo.workflow import BuildVG, Findmotif
from grafimo.constructVG import construct_vg, indexVG
from grafimo.motif import get_motif_pwm
from grafimo.motif_set import MotifSet
from grafimo.extract_regions import get_regions
from grafimo.res_writer import print_results, write_results
from grafimo.score_sequences import compute_results
import time
import pandas as pd

# version of GRAFIMO
__version__ = '1.0.1'


def buildvg(args_obj):

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


def findmotif(args_obj):

    if not isinstance(args_obj, Findmotif):
        raise ValueError("Unknown arguments object type")

    printWelcomeMsg()

    # if verbose == True will be printed a lot of unusefull stuff
    verbose = args_obj.get_verbose()

    cores = args_obj.get_cores()  # cores to use

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
        print("\t- Top graphs: ", args_obj.get_top_graphs())
        print("\t- no-qvalue: ", args_obj.get_no_qvalue())
        print("\t- no-reverse: ", args_obj.get_no_reverse())
        print("\t- Text only: ", args_obj.get_text_only())
        print("\t- Threshold on q-values: ", args_obj.get_qvalueT())
        print("\t- Verbose: ", verbose)
        print("\t- Test mode: ", args_obj.get_test())
        print()  # newline
    # end if

    # if we have a whole genome graph, check that it is indexed (XG)
    # index it, otherwise
    if args_obj.has_graph_genome():
        vg = args_obj.get_graph_genome()

        if vg.split('.')[-1] == 'vg':
            code = indexVG(vg)

            if code != 0:
                errmsg = ''.join(["\n\nERROR: unable to index ", vg, ". Exiting"])
                raise(errmsg)
            # end if

            vg = vg.split('.vg')[-2]
            vg = ''.join([vg, '.xg'])
            args_obj.set_xg(vg)
        # end if
    # end if

    # process the motifs
    motifs = args_obj.get_motif()
    mtfSet = MotifSet()

    for mtf in motifs:
        if verbose:
            print("Processing motif(s) in: ", mtf, "\n")

        m = get_motif_pwm(mtf, args_obj, cores)

        if not args_obj.get_test():
            for i in range(len(m)):
                print("Score matrix:")
                m[i].print("score_matrix")
        # end if

        mtfSet.addMotif(m)
    # end for

    motif_num = mtfSet.length()

    for mtf in mtfSet:
        # extract sequences
        sequence_loc = get_regions(mtf, args_obj)

        # score sequences
        res = compute_results(mtf, sequence_loc, args_obj)

        # write results
        if args_obj.get_text_only():
            print_results(res)
        else:
            write_results(res, mtf.getMotifID(), motif_num, args_obj)
        # end if
    # end for

# end of findmotif()


def printWelcomeMsg():
    """
        Prints the intro message for GRAFIMO
        ----
        Params:
            None
        ----
        Returns:
            None
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

