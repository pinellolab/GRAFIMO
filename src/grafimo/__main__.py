#!/usr/bin/env python
#
# MIT License
#
# Copyright (c) 2020 Pinello Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


"""

GRAFIMO version {version}

Copyright (C) 2020 Manuel Tognon <manu.tognon@gmail.com> <manuel.tognon@studenti.univr.it>

GRAph-based Find of Individual Motif Occurrences.

GRAFIMO scans a given genome variation graphs for occurrences of a given motif, searching in
a given set opf genomic coordinates, defined as a BED file.
GRAFIMO allows the user, also, to build from scratch a set of genome variation graphs for
each chromosome or of a user defined subset of them.

Usage:

    - to scan a precomputed genome variation graph:
        grafimo findmotif -g GENOME_GRAPH.xg -b BEDFILE.bed -m MOTIF.meme [options]

    - to scan a set of precomputed genome variation graphs:
        grafimo findmotif -d path/to/my/vgs -b BEDFILE.bed -m MOTIF.meme [options]

    - to build a new set of genome variation graphs (one for each chromosome or a user
      defined subset):
        grafimo buildvg -l REFERENCE.fa -v VCF.vcf.gz [options]

GRAFIMO will scan the given genome variation graph or all the ones contained in a given
directory, to search the given motif occurrences in the regions defined by the user in
the given BED file. Once GRAFIMO finishes to scan the VG (or VGs) it stores the results
in an output directory (by default grafimo_out_JOBID_MOTIFID) in three different file
format:
    - TSV
    - HTML
    - GFF3

See https://github.com/pinellolab/GRAFIMO or https://github.com/InfOmics/GRAFIMO for the full documentation.

If you use GRAFIMO in your research, please cite us:

    Tognon, M., et al., "GRAFIMO: variant-aware motif scanning via genome variation graphs", xxxx


Run 'grafimo --help'to see all command-line options.
"""


from argparse import ArgumentParser, SUPPRESS, HelpFormatter
from grafimo.workflow import BuildVG, Findmotif
from grafimo.utils import die, initialize_chroms_list, isJaspar_ff, isMEME_ff, check_deps, sigint_handler, \
        EXT_DEPS, CHROMS_LIST
from grafimo.grafimo import __version__, buildvg, findmotif
from grafimo.GRAFIMOException import DependencyError
import multiprocessing as mp
import glob
import time
import sys
import os


class GRAFIMOArgumentParser(ArgumentParser):
    """

        This is a redefinition of the arguments parser:
            - The usage message is not prefixed by 'usage'
            - A brief message is printed instead of the full usage message

    """

    class GRAFIMOHelpFormatter(HelpFormatter):

        def add_usage(self, usage, actions, groups, prefix='None'):
            if usage is not SUPPRESS:
                args = usage, actions, groups, ''
                self._add_item(self._format_usage, args)

    # end of GRAFIMOHelpFormatter

    def __init__(self, *args, **kwargs):
        kwargs['formatter_class'] = self.GRAFIMOHelpFormatter
        kwargs['usage'] = kwargs['usage'].replace("{version}", __version__)
        super().__init__(*args, **kwargs)

    def error(self, msg):
        sys.stderr.write("Run 'grafimo --help' to see the usage\n")
        self.exit(2, "\n{0}: ERROR: {1}\n".format(self.prog, msg))

# end of GRAFIMOArgumentParser


class CmdError(Exception):
    pass


def get_AP():
    """
        Get the arguments from the command-line call of
        GRAFIMO and return the parser containing all the input arguments
        ----
        Parameters:
            None
        ----
        Return:
            parser (GRAFIMOArgumentParser)
    """

    # define the argument parser
    parser = GRAFIMOArgumentParser(usage=__doc__, add_help = False)

    #####################################################################
    # start reading the arguments
    #####################################################################

    group = parser.add_argument_group("Options")
    group.add_argument("-h", "--help", action="help", help="Show the help message and exit.")
    group.add_argument("--version", action="version", help="Show the current version and exit.", version=__version__)
    group.add_argument("--cores", type=int, default=0, metavar="NCORES", nargs='?', const=0,
                       help="Number of cores to use. The default value is %(default)s (auto-detection)."
                            "If you chose to query the whole genome variation graph note that the default "
                            "option is to use only one core to avoid memory issues; if you are running "
                            "GRAFIMO on a machine with RAM >= 128 GB maybe you can use more cores.")

    #####################################################################
    # mandatory arguments
    #####################################################################
    group.add_argument("workflow", type=str, default='',
                            help="This option accpets only two values: 'buildvg' or 'findmotif'. "
                                 "This option is mandatory. When used 'buildvg', GRAFIMO will "
                                 "compute the genome variation graph for each chromosome or a user "
                                 "defined subset, from the given reference genome and VCF file. "
                                 "When used 'findmotif' GRAFIMO will scan a given VG or all the VGs "
                                 "contained in a given directory, for the given motif occurrences")

    # arguments to build the VGs
    group.add_argument("-l", "--linear-genome", type=str, help='Path to the linear genome (FASTA format required)',
                       nargs='?', default='', metavar='LINEAR-GENOME', dest='linear_genome')
    group.add_argument("-v", "--vcf", type=str, nargs='?', default='', metavar='VCF',
                       help='Path to the VCF file. The vcf must be compressed '
                            '(e.g. myvcf.vcf.gz)')

    # arguments to scan the VG/VGs
    group.add_argument("-g", "--graph-genome", type=str, nargs='?', metavar='GRAPH-GENOME',
                       help='Path to the VG genome graph file (VG or XG format)',
                       default='', dest='graph_genome')
    group.add_argument("-d", "--graph-genome-dir", type=str, nargs='?', metavar='GRAPH-GENOMES-DIR',
                       help='Path to a directory containing a variable number of VGs '
                            'graph genomes (VG or XG format)',
                       default='', dest='graph_genome_dir')
    group.add_argument("-b", "--bedfile", type=str, help='Path to the BED file containing the regions to scan',
                       metavar='BEDFILE')
    group.add_argument("-m", "--motif", type=str, nargs='+', metavar=('MOTIF1', 'MOTIF2'),
                       help='Path to the motif file (JASPAR or MEME format required)')

    #####################################################################
    # optional arguments
    #####################################################################
    group.add_argument("-c", "--chroms", type=str, nargs='*',
                       help='Chromosomes for which the VG will be built or in which GRAFIMO '
                            'will search the given motif.\n'
                            'To consider all the chromosomes, just skip this argument.\n'
                            'This argument can be used during the building of VGs or during '
                            'their scan for the occurrences of the given motif',
                       default=[], metavar=('1', 'X'))
    group.add_argument("-k", "--bgfile", type=str,
                       help='Path to the background file [optional]', nargs='?',
                       const='', default='UNIF', metavar='BACKGROUND')
    group.add_argument("-p", "--pseudo", type=float,
                       help='Pseudocount to add to motif counts [optional]', nargs='?', const='0.1',
                       default='0.1', metavar='PSEUDOCOUNT')
    group.add_argument("-t", "--threshold", type=float, nargs='?', default=1e-4,
                       metavar='THRESHOLD', const='1e-4',
                       help='Threshold that will be applied on the P-values (by default) of each found '
                            'motif occurrence. It is possible to apply the threshold on the q-values '
                            'using the -q (--qvalueT) option.\n'
                            'Default is 1e-4 [optional]')
    group.add_argument("-q", "--no-qvalue", action='store_true', default=False, dest='no_qvalue',
                       help='With this parameter the q-values will not be computed')
    group.add_argument("-r", "--no-reverse", default=False, action='store_true', dest='no_reverse',
                       help='With this parameter GRAFIMO will scan only the forward strand')
    group.add_argument("-f", "--text-only", default=False, action='store_true',
                       dest='text_only',
                       help="Print the results in TSV directly to the standard output",
                       )
    group.add_argument("-o", "--out", type=str,
                       help='Name of the directory where the results will be stored [optional]',
                       nargs='?', const='grafimo_out', default='grafimo_out', metavar='OUTDIR')
    group.add_argument("--qvalueT", action='store_true', default=False, dest='qval_t',
                       help='The threshold will be applied on the q-values, rather than the P-values')
    group.add_argument("--top-graphs", type=int,
                       help='The PNG image of the regions containing the top GRAPHS_NUM sequences '
                            '(sorted by P-value) will be stored in the output directory',
                       nargs='?', const=0, default=0, metavar='GRAPHS_NUM', dest='top_graphs')

    group.add_argument("--verbose", default=False, action='store_true',
                       help="Output a lot of additional information about the execution")

    return parser

# end of get_AP()


def main(cmdLineargs=None):
    """

        Main function of GRAFIMO.

        The arguments given in input are checked for consistency,
        then a pipeline is followed.

        ----
        Parameters:
            cmdLineargs (str)
        ----
        Returns:
            None

    """

    try:
        # starting point of the execution time
        start = time.time()

        # read the command-line arguments
        parser = get_AP()

        if cmdLineargs is None:
            cmdLineargs = sys.argv[1:] # take input args

        # the second argument must be buildvg or findmotif
        if ((cmdLineargs[0] != "-h" and cmdLineargs[0] != "--help") and
                cmdLineargs[0] != "--version" and
                (cmdLineargs[0] != "buildvg" and cmdLineargs[0] != "findmotif")):
            parser.error("The second argument must be one between 'buildvg' and 'findmotif'")
            die(1)

        args = parser.parse_args(cmdLineargs)  # parse args

        if args.verbose:
            print("Parsing arguments...")
            start_args_parse = time.time()

        #####################################################################
        # check arguments consistency
        #####################################################################

        if args.workflow != "buildvg" and args.workflow != "findmotif":
            parser.error("Do not know what to do. Available options: create VGs with 'grafimo buildvg' or scan a "
                         "precomputed genome variation graph with 'grafimo findmotif'")
            die(1)

        # cores (shared by the two workflows)
        if args.cores < 0:
            parser.error("The number of cores cannot be negative")

        elif args.cores == 0 and args.graph_genome:
            args.cores = 1     # to query a whole genome graph is loaded into RAM, since usually are
                                 # very heavy in terms of bytes is safer to use 1 thread by default, otherwise
                                 # it would be loaded #cores times. If you want use more cores, be sure
                                 # your system can handle the resulting amount of data

        elif args.cores == 0:
            args.cores = mp.cpu_count()  # by default take all the available CPUs
        # end if

        # check verbose flag
        if (not isinstance(args.verbose, bool) or
                (args.verbose != False and args.verbose != True)):
            parser.error('The --verbose parameter accepts only True or False values')

        # chromosomes check (shared by the two workflows)
        for c in args.chroms:
            if c not in CHROMS_LIST:
                parser.error("Invalid chromosome")
                
        args.chroms = initialize_chroms_list(args.chroms)

        # checks for buildvg workflow
        if args.workflow == "buildvg":

            if args.graph_genome_dir:
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif args.graph_genome:
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif args.bedfile:
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif args.motif:
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif args.bgfile != 'UNIF':  # if default ignored
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif args.pseudo != 0.1:  # if default ignored
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif args.threshold != 1e-4:  # if default ignored"
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif args.no_qvalue:
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif args.no_reverse:
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif args.text_only:
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif args.qval_t:
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif args.top_graphs != 0:  # if default ignored
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif not args.linear_genome:
                parser.error("No reference genome given")
                die(1)

            elif not args.vcf:
                parser.error("No VCF file given")
                die(1)

            else:
                # check linear genome
                if (args.linear_genome.split('.')[-1] != 'fa' and
                        args.linear_genome.split('.')[-1] != 'fasta'):
                    parser.error('The linear genome must be in FASTA format (FASTA and FA extensions allowed)')
                    die(1)

                else:
                    if len(glob.glob(args.linear_genome)) != 1:
                        parser.error('Cannot find the given reference genome file')
                        die(1)

                    args.linear_genome = os.path.abspath(args.linear_genome)
                # end if

                # check VCF
                if ((args.vcf.split('.')[-1] != 'gz' and args.vcf.split('.')[-1] != 'zip')
                        or args.vcf.split('.')[-2] != 'vcf'):  # allow only compressed VCF files
                    parser.error('Incorrect VCF file given: the VCF must be compressed (e.g. myvcf.vcf.gz)')
                    die(1)

                else:
                    if len(glob.glob(args.vcf)) <= 0:
                        parser.error('Cannot find the given VCF file')
                        die(1)
                    args.vcf = os.path.abspath(args.vcf)

                # by deafult the built VGs will be stored in the current directory
                if args.out == "grafimo_out":  # general default value
                    args.out = os.path.abspath("./")

                workflow = BuildVG(args)

                if args.verbose:
                    end_args_parse = time.time()
                    print(''.join(["Arguments parsed in ", str(end_args_parse - start_args_parse), "s"]))
            # end if
        # end if

        # checks for findmotif workflow
        if args.workflow == "findmotif":
            if args.linear_genome:
                parser.error("Invalid arguments for grafimo findmotif")
                die(1)

            elif args.vcf:
                parser.error("Invalid arguments for grafimo buildvg")
                die(1)

            elif not args.graph_genome_dir and not args.graph_genome:
                parser.error("No genome variation graph or directory containing them given")
                die(1)

            elif not args.bedfile:
                parser.error("No BED file given")
                die(1)

            elif not args.motif:
                parser.error("No motif file (MEME of JASPAR format) given")
                die(1)

            else:

                # only one between graph_genome and graph_genome_dir allowed
                if args.graph_genome and args.graph_genome_dir:
                    parser.error("Invalid arguments for grafimo buildvg")
                    die(1)

                # check graph_genome
                if args.graph_genome:
                    if (args.graph_genome.split('.')[-1] != 'xg' and
                            args.graph_genome.split('.')[-1] != 'vg'):
                        parser.error("Cannot use the given genome variation graph (only VG or XG format allowed)")
                        die(1)

                    elif not os.path.isfile(args.graph_genome):
                        parser.error("Unable to find the given variation genome graph")
                        die(1)

                    else:
                        graph_genome = os.path.abspath(args.graph_genome)  # safer to use absolute path
                        args.graph_genome = graph_genome
                    # end if
                # end if

                # check graph_genome_dir
                if args.graph_genome_dir:
                    if not os.path.isdir(args.graph_genome_dir):
                        parser.error("Cannot find the given directory containing the genome variation graphs")
                        die(1)

                    if args.graph_genome_dir[-1] == '/':
                        graph_genome_dir = args.graph_genome_dir

                    else:
                        graph_genome_dir = ''.join([args.graph_genome_dir, '/'])
                    # end if

                    if len(glob.glob(graph_genome_dir + '*.xg')) <= 0:
                        parser.error(' '.join(['No XG genome variation graph found in', graph_genome_dir]))
                        die(1)

                    else:
                        graph_genome_dir = os.path.abspath(graph_genome_dir)
                        args.graph_genome_dir = graph_genome_dir
                    # end if
                # end if

                # check BED file
                if args.bedfile:
                    if args.bedfile.split('.')[-1] != 'bed':
                        parser.error('Incorrect BED file given')
                        die(1)

                    else:
                        bedfile = args.bedfile

                        if len(glob.glob(bedfile)) <= 0:
                            parser.error('Cannot find the given BED file')
                    # end if

                else:
                    parser.error('No BED file given')
                # end if

                # check motif file
                if not args.motif:
                    parser.error('No motif given')

                else:
                    motifs = args.motif

                    # check if the given motifs exist
                    for m in motifs:
                        if not isMEME_ff(m) and not isJaspar_ff(m):
                            parser.error("Unrecognized motif file format (only MEME or JASPAR allowed)")
                            die(1)

                        if len(glob.glob(m)) <= 0:
                            parser.error('Cannot find motif file: ' + m)
                            die(1)
                    # end for
                # end if

                # check background file
                if args.bgfile != 'UNIF':
                    bgfile = args.bgfile  # we have a path to a bg file

                    if len(glob.glob(bgfile)) <= 0:
                        parser.error('Cannot find the given background file')
                        die(1)
                # end if

                # check pseudocount
                if args.pseudo <= 0:
                    parser.error('The pseudocount cannot be less than or equal 0')
                    die(1)

                # check threshold
                if args.threshold <= 0 or args.threshold > 1:
                    parser.error('The pvalue threshold must be between 0 and 1')
                    die(1)

                # check q-value flag
                if (not isinstance(args.no_qvalue, bool) or
                        (args.no_qvalue != False and args.no_qvalue != True)):
                    parser.error('The --qvalue parameter accepts only True or False as values')
                    die(1)

                # check no reverse flag
                if (not isinstance(args.no_reverse, bool) or
                        (args.no_reverse != False and args.no_reverse != True)):
                    parser.error('The --no-reverse parameter accepts only True or False as values')
                    die(1)

                # check text only flag
                if (not isinstance(args.text_only, bool) or
                        (args.text_only != False and args.text_only != True)):
                    parser.error('The --text-only parameter accepts only True or False values')
                    die(1)

                # out directory
                if args.out == 'grafimo_out':  # default option
                    # to make unique the output directory we add the PID
                    # to the name.
                    #
                    # This is useful when calling grafimo in different runs on the
                    # same machine.

                    args.out = ''.join([args.out, '_', str(os.getpid())])

                # check threshold on q-value flag
                if (not isinstance(args.qval_t, bool) or
                        (args.qval_t != False and args.qval_t != True)):
                    parser.error("The --qvalueT parameter accepts only True or False as values")
                    die(1)

                elif args.no_qvalue == True and args.qval_t == True:
                    parser.error("Cannot apply the threshold on q-values if you don't want them")
                    die(1)

                # check the number of graph regions to store as PNG images
                if args.top_graphs < 0:
                    parser.error("The number of region graphs to show must be positive")

                workflow = Findmotif(args)

                if args.verbose:
                    end_args_parse = time.time()
                    print(''.join(["Arguments parsed in ", str(end_args_parse - start_args_parse), "s"]))

            # end if
        # end if

        # check that external dependencies are satisfied
        if args.verbose:
            print("Checking GRAFIMO external dependencies " + str(EXT_DEPS))
            start_deps = time.time()

        satisfied, deps_lack = check_deps()

        if not satisfied and len(deps_lack) > 0:
            raise DependencyError("\n\nERROR: The following dependencies are not sastisfied: " +
                                      str(deps_lack) +
                                      "\nPlease, solve them before running GRAFIMO")
            die(1)

        elif not satisfied and len(deps_lack) <= 0:
            raise DependencyError("Some dependencies were found, but was not possible to track them." 
                                        "\nBe sure they are available in system PATH")
            die(1)
        # end if

        if args.verbose and satisfied:
            end_deps = time.time()
            print("Dependencies correctly satisfied")
            print(''.join(["Dependencies checked in ", str(end_deps - start_deps), "s"]))

        #####################################################################

        """
            dependency check was ok, so we go to workflow selection:
               - creation of the genome variation graph for 
                   each chromosome or a user defined subset of them
               - scan of a precomputed VG or a set of precomputed VG
        """

        if isinstance(workflow, BuildVG):
            # build the VG for each chromosome or a user defined subset of them
            buildvg(workflow)

        elif isinstance(workflow, Findmotif):
            # scan a precomputed VG or a set of VGs
            findmotif(workflow)

        else:
            raise ValueError("Unknown arguments object type")
        # end if

        end = time.time()  # GRAFIMO execution finishes here

        print(''.join(["\nElapsed time: ", str(end - start), "s"]))

    except KeyboardInterrupt:
        sigint_handler()

    finally:
        pass
    # end try

# end of main()


#####################################################################
#
# Entry point for GRAFIMO
#
#####################################################################

if __name__ == "__main__":
    main()

