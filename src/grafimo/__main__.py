#!/usr/bin/env python
#
# MIT License
#
# Copyright (c) 2021 Pinello Lab
#
# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the 
# "Software"), to deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, merge, publish, 
# distribute, sublicense, and/or sell copies of the Software, and to 
# permit persons to whom the Software is furnished to do so, subject to 
# the following conditions:
#
# The above copyright notice and this permission notice shall be 
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


"""
GRAFIMO version {version}

Copyright (C) 2021 Manuel Tognon <manu.tognon@gmail.com> <manuel.tognon@univr.it>

GRAph-based Find of Individual Motif Occurrences.

GRAFIMO scans genome variation graphs in VG format (Garrison et al., 2018) 
to find potential binding site occurrences of DNA motif(s).

Usage:

    * to scan a set of precomputed genome variation graphs (-- SUGGESTED --):

        grafimo findmotif -d path/to/my/vgs -b BEDFILE.bed -m MOTIF.meme [options]

    * to scan a precomputed genome variation graph:
        
        grafimo findmotif -g GENOME_GRAPH.xg -b BEDFILE.bed -m MOTIF.meme [options]

    * to build a new set of genome variation graphs (one for each chromosome or a user defined subset):
        
        grafimo buildvg -l REFERENCE.fa -v VCF.vcf.gz [options]

GRAFIMO searches DNA motif(s) occurrences, given as a PWM (MEME, JASPAR format), 
in a set of genomic coordinates on VGs. Genomic coordinates are provided in tab-separated file
(e.g. ENCODE narrowPeak or UCSC BED format).

GRAFIMO scans in a single run all the genomes encoded in the input VG. 
GRAFIMO helps studying how genetic variation present within a 
population of individuals affects the binding affinity score of input DNA 
motif(s).

GRAFIMO provides the possibility to build from user data a new genome 
variation graph, by building a VG for each chromosome.

GRAFIMO results are reported in three files (stored in a directory named `grafimo_out_JOBID` by default):
    
    * TSV report, containing all the reported binding site candidates, with 
      an associated score, P-value, q-value, haplotype frequency within 
      the analyzed population and a flag stating if the current occurrence 
      contains a genomic variants

    * HTML version of the TSV report

    * GFF3 file to load the found motif candidates on the UCSC genome browser 
      as custom track

Citation:

    Tognon, Manuel, et al. "GRAFIMO: variant and haplotype aware motif scanning on 
    pangenome graphs." bioRxiv (2021).


Run "grafimo --help" to see all command-line options.
See https://github.com/pinellolab/GRAFIMO/wiki or https://github.com/InfOmics/GRAFIMO/wiki for the full documentation.
"""


from grafimo.GRAFIMOArgumentParser import GRAFIMOArgumentParser
from grafimo.workflow import BuildVG, Findmotif
from grafimo.utils import die, initialize_chroms_list, isJaspar_ff, isMEME_ff, \
    check_deps, sigint_handler, EXT_DEPS, CHROMS_LIST, DEFAULT_OUTDIR, NOMAP, ALL_CHROMS, anydup, exception_handler, UNIF, isbed
from grafimo.grafimo import __version__, buildvg, findmotif
from grafimo.GRAFIMOException import DependencyError

from typing import List, Optional
from glob import glob

import multiprocessing as mp

import argparse
import time
import sys
import os


def get_parser() -> GRAFIMOArgumentParser:
    """Read the comman-line arguments given to GRAFIMO and creates a 
    GRAFIMOArgumentParser object containing all of them

    Parameters
    ----------

    Returns
    -------
    GRAFIMOArgumentParser
        parser object containing the command-line arguments given to 
        GRAFIMO
    """

    # define the argument parser
    parser: GRAFIMOArgumentParser = GRAFIMOArgumentParser(usage=__doc__, 
                                                          add_help = False)

    # start commandline arguments parsing

    # general options
    group = parser.add_argument_group("General options")
    group.add_argument("-h", "--help", action="help", help="Show this help " 
                       "message and exit.")
    group.add_argument("--version", action="version", help="Show software "
                       "version and exit", version=__version__)
    group.add_argument("-j", "--cores", type=int, default=0, metavar="NCORES", 
                       nargs='?', const=0, help="Number of CPU cores to use. Use 0 "
                       "to auto-detect. Default: %(default)s. "
                       "To search motifs in a whole genome variation graph the "
                       "default is 1 (avoid memory issues).")
    group.add_argument("--verbose", default=False, action='store_true',
                       help="Print additional information about GRAFIMO run.")
    group.add_argument("--debug", action="store_true", default=False, dest="debug",
                       help="Enable error traceback.")
    group.add_argument("-o", "--out", type=str, help="Output directory."
                       " [optional]", nargs='?', const='grafimo_out', default="", 
                       metavar='OUTDIR')
    group.add_argument("workflow", type=str, default='',
                       help="Mandatory argument placed immediately after \"grafimo\". "
                       "Only two values are accepted: \"buildvg\" and \"findmotif\". "
                       "When called \"grafimo buildvg\", the software will compute "
                       "the genome variation graph from input data "
                       "(see \"buildvg options\" section below for more arguments). "
                       "When called \"grafimo findmotif\", the software will scan "
                       "the input VG(s) for potential occurrences of the input "
                       "motif(s) (see \"findmotif option\" section below for more arguments)."
                       )

    # buildvg options
    group = parser.add_argument_group("Buildvg options")
    group.add_argument("-l", "--linear-genome", type=str, nargs='?', 
                       default='', metavar='REFERENCE-FASTA', dest='linear_genome',
                       help="Path to reference genome FASTA file.")
    group.add_argument("-v", "--vcf", type=str, nargs='?', default='', 
                       metavar='VCF', help="Path to VCF file. Note that the VCF "
                       "should be compressed (e.g. myvcf.vcf.gz).")
    group.add_argument("--chroms-build", type=str, nargs="*", default=[], 
                        metavar=("1", "X"), dest="chroms_build",
                        help="Chromosomes for which construct the VG. By default "
                        "GRAFIMO constructs the VG for all chromsomes.")
    group.add_argument("--chroms-prefix-build", type=str, nargs="?", default="",
                       metavar="CHRPREFIX", dest="chroms_prefix_build",
                       help="Prefix to append in front of chromosome numbers. "
                       "To name chromosome VGs with only their number (e.g. 1.xg), "
                       "use \"--chroms-prefix-build \"\" \". Default: %(default)s.")
    group.add_argument("--chroms-namemap-build", type=str, nargs="?", default="NOMAP",
                      metavar="NAME-MAP-FILE", dest="chroms_namemap_build",
                      help="Space or tab-separated file, containing original "
                      "chromosome names in the first columns and the names to "
                      "use when storing corresponding VGs. By default the VGs "
                      "are named after the encoded chromosome (e.g. chr1.xg).")
    group.add_argument("--reindex", action='store_true', default=False, 
                       help="Reindex the VCF file with Tabix, even if a TBI "
                       "index os already available.")

    # findmotif options
    group = parser.add_argument_group("Findmotif options")
    group.add_argument("-g", "--genome-graph", type=str, nargs='?', 
                       metavar='GENOME-GRAPH', default='', dest='graph_genome',
                       help="Path to VG pangenome variation graph (VG or XG "
                       "format).")
    group.add_argument("-d", "--genome-graph-dir", type=str, nargs='?', 
                       metavar='GENOME-GRAPHS-DIR', default='', 
                       dest='graph_genome_dir', help="Path to the directory "
                       "containing the pangenome variation graphs to scan (VG "
                       "or XG format)")
    group.add_argument("-b", "--bedfile", type=str, metavar='BEDFILE',
                        help="BED file containing the genomic regions to scan "
                        "for occurrences of the input motif(s).")
    group.add_argument("-m", "--motif", type=str, nargs='+', metavar=('MOTIF1', 
                       'MOTIF2'), help="Motif Position Weight Matrix (MEME or "
                       "JASPAR format).")
    group.add_argument("-k", "--bgfile", type=str, nargs='?', const='', 
                       default=UNIF, metavar='BACKGROUND', help="Background "
                       "distribution file.")
    group.add_argument("-p", "--pseudo", type=float, nargs='?', const='0.1',
                       default='0.1', metavar='PSEUDOCOUNT', help="Pseudocount "
                       "used during motif PWM processing.")
    group.add_argument("-t", "--threshold", type=float, nargs='?', default=1e-4,
                       metavar='THRESHOLD', const='1e-4', help="Statistical "
                       "significance threshold value. By default the threshold "
                       "is applied on P-values. To apply the threshold on "
                       "q-values use the \"--qvalueT\" options. Default:"
                       "%(default)s.")
    group.add_argument("-q", "--no-qvalue", action='store_true', default=False, 
                       dest='no_qvalue', help="If used, GRAFIMO skips q-value "
                       "computation.")
    group.add_argument("-r", "--no-reverse", default=False, action='store_true', 
                       dest='no_reverse', help="If used, GRAFIMO scans only the "
                       "forward strand.")
    group.add_argument("-f", "--text-only", default=False, action='store_true',
                       dest='text_only', help="Print results to stdout.")
    group.add_argument("--chroms-find", type=str, nargs="*", default=[],
                       metavar=("1", "X"), dest="chroms_find", help="Scan only "
                       "the specified chromosomes.")
    group.add_argument("--chroms-prefix-find", type=str, nargs="?", default="",
                       metavar="CHRPREFIX", dest="chroms_prefix_find", 
                       help="Prefix shared by all chromosomes. The prefix should "
                       "be followed by the chromosome number. If chromosome VGs "
                       "are stored only with their chromosome number (e.g. 1.xg) "
                       "use \"--chroms-prefix-fin \"\" \". Default: %(default)s.")
    group.add_argument("--chroms-namemap-find", type=str, nargs="?", default="NOMAP",
                       metavar="NAME-MAP-FILE", dest="chroms_namemap_find",
                       help="Space or tab-separated file, containing original "
                      "chromosome names in the first columns and the names used "
                      "to store the corresponding VGs. By default GRAFIMO assumes "
                      "that VGs are named after the encoded chromosome (e.g. chr1.xg).")
    group.add_argument("--recomb", action='store_true', default=False, 
                       help="Consider all the possible recombinants sequences "
                       "which could be obtained from the genetic variants encoded "
                       "in the VG. With this option the haplotypes encoded in "
                       "the VG are ignored.")
    group.add_argument("--qvalueT", action='store_true', default=False, 
                       dest='qval_t', help="Apply motif occurrence score "
                       "statistical significance threshold on q-values rather "
                       "than on P-values.")
    group.add_argument("--top-graphs", type=int, nargs='?', const=0, default=0, 
                       metavar='GRAPHS-NUM', dest='top_graphs', help="Store the "
                       "PNG image of the top \"GRAPHS-NUM\" regions of the VG "
                       "(motif occurrences sorted by increasing P-value).")

    return parser

# end of get_parser()


def main(cmdLineargs: Optional[List[str]] = None) -> None :

    try:
        # starting point of the execution time
        start: float = time.time()

        # read the command-line arguments
        parser: GRAFIMOArgumentParser = get_parser()

        if cmdLineargs is None:
            cmdLineargs: List[str] = sys.argv[1:]  # get input args

        # no arguments given --> print help
        if len(cmdLineargs) == 0:
            parser.error_noargs()
            die(2)

        # the second argument must be buildvg or findmotif
        if ((cmdLineargs[0] != "-h" and cmdLineargs[0] != "--help") and
                cmdLineargs[0] != "--version" and
                (cmdLineargs[0] != "buildvg" and cmdLineargs[0] != "findmotif")):
            parser.error(
                "The second argument must be one between 'buildvg' and 'findmotif'")
            die(1)

        args: argparse.Namespace = parser.parse_args(cmdLineargs) 

        if args.verbose:
            print("Parsing arguments...")
            start_args_parse: float = time.time()

        #--------------------------------------------------------------#
        # check commandline arguments consistency
        #

        #---------------------- general options -----------------------#
        
        # workflow type
        if args.workflow != "buildvg" and args.workflow != "findmotif":
            parser.error(
                "Unexpected workflow given. Available options:\n"
                "\tbuildvg: construct VG from user data.\n"
                "\tfindmotif: scan VG for DNA motif(s) occurrences")
            die(1)

        # cpu cores 
        if args.cores < 0:
            parser.error("Negative number of CPU cores given")
        elif args.cores == 0 and args.graph_genome:
            # when whole genome variation graph is given, it is safer to
            # use 1 CPU core by default. This beacuse of the space needed
            # to load the whole VG on RAM.
            #
            # CAVEAT: before requiring more CPU cores to be used, be sure
            # your system has enough memory  
            args.cores = 1  
        elif args.cores == 0:
            # default option -> use all available CPU cores
            args.cores = mp.cpu_count() 
        else:  # args.cores > 0
            if args.cores > mp.cpu_count():
                parser.error(
                    "Too many CPU cores to use ({})".format(args.cores)
                )

        # verbosity 
        if (not isinstance(args.verbose, bool) or
                (args.verbose != False and args.verbose != True)):
            parser.error(
                '\"--verbose\" does not accept any positional argument'
            )

        # debugging
        if (not isinstance(args.debug, bool) or
                (args.debug != False and args.debug != True)):
            parser.error(
                "\"--debug\" does not accept any positional argument"
            )

        #---------------------- buildvg options -----------------------#

        buildvg_err_msg: str = "Unexpected arguments for \"grafimo buildvg\": \"{}\""

        if args.workflow == "buildvg":

            if args.graph_genome_dir:
                parser.error(buildvg_err_msg.format("-d, --genome-graph-dir"))
                die(1)
            elif args.graph_genome:
                parser.error(buildvg_err_msg.format("-g, --genome-graph"))
                die(1)
            elif args.bedfile:
                parser.error(buildvg_err_msg.format("-b, --bedfile"))
                die(1)
            elif args.motif:
                parser.error(buildvg_err_msg.format("-m, --motif"))
                die(1)
            elif args.bgfile != UNIF:  # if default ignored
                parser.error(buildvg_err_msg.format("-k, --bgfile"))
                die(1)
            elif args.pseudo != 0.1:  # if default ignored
                parser.error(buildvg_err_msg.format("-p, --pseudo"))
                die(1)
            elif args.threshold != 1e-4:  # if default ignored
                parser.error(buildvg_err_msg.format("-t, --thresh"))
                die(1)
            elif args.no_qvalue:
                parser.error(buildvg_err_msg.format("-q, --no-qvalue"))
                die(1)
            elif args.no_reverse:
                parser.error(buildvg_err_msg.format("-r, --no-reverse"))
                die(1)
            elif args.text_only:
                parser.error(buildvg_err_msg.format("-f, --text-only"))
                die(1)
            elif args.chroms_find:
                parser.error(buildvg_err_msg.format("--chroms-find"))
                die(1)
            elif args.chroms_prefix_find: 
                parser.error(buildvg_err_msg.format("--chroms-prefix-find"))
                die(1)
            elif args.chroms_namemap_find != NOMAP:  # if default ignored
                parser.error(buildvg_err_msg.format("--chroms-namemap-find"))
                die(1)
            elif args.qval_t:
                parser.error(buildvg_err_msg.format("--qvalueT"))
                die(1)
            elif args.recomb:
                parser.error(buildvg_err_msg.format("--recomb"))
                die(1)
            elif args.top_graphs != 0:  # if default ignored
                parser.error(buildvg_err_msg.format("--top-graphs"))
                die(1)
            elif not args.linear_genome:
                parser.error("No reference genome given")
                die(1)
            elif not args.vcf:
                parser.error("No VCF file given")
                die(1)
            else:  # arguments for buildvg are correct
                # reference genome
                if (args.linear_genome.split('.')[-1] != 'fa' and
                        args.linear_genome.split('.')[-1] != 'fasta'):
                    parser.error(
                        "The reference genome file must be in FASTA format"
                    )
                    die(1)
                else:
                    if not os.path.isfile(args.linear_genome):
                        parser.error(
                            "Unable to find {}".format(args.linear_genome)
                        )
                        die(1)
                    if os.stat(args.linear_genome).st_size == 0:  # empty file
                        parser.error(
                            "{} seems to be empty.".format(args.linear_genome)
                        )
                        die(1)
                    args.linear_genome = os.path.abspath(args.linear_genome)
                # VCF --> the VCF file must have been compressed with
                # bgzip (https://github.com/samtools/tabix)
                if (args.vcf.split(".")[-1] != "gz" and 
                    args.vcf.split(".")[-2] != "vcf"):  
                    parser.error(
                        "Wrong VCF file given. The VCF file must have been "
                        "compressed with bgzip (e.g. myvcf.vcf.gz)"
                    )
                    die(1)
                else:
                    if not os.path.isfile(args.vcf):
                        parser.error('Unable to find {}'.format(args.vcf))
                        die(1)
                    if os.stat(args.vcf).st_size == 0:  # empty file
                        parser.error(
                            "{} seems to be empty.".format(args.vcf)
                        )
                        die(1)
                    args.vcf = os.path.abspath(args.vcf)

                # chromosome to construct VG
                if len(args.chroms_build) == 0:
                    args.chroms_build = [ALL_CHROMS]  # use all chromosome 
                else:
                    if anydup(args.chroms_build):
                        parser.error(
                            "Duplicated chromosome names given to \"--chroms-build\""
                        )

                # chromosome name-map
                if args.chroms_namemap_build != NOMAP:
                    if not os.path.isfile(args.chroms_namemap_build):
                        parser.error(
                            "Unable to locate {}".format(args.chroms_namemap_build)
                        )
                if (args.chroms_prefix_build and args.chroms_namemap_build != NOMAP):
                    parser.error(
                        "\"--chroms-prefix-build\" and \"chroms-namemap-build\" "
                        "cannot used together. Choose one of those options"
                    )

                # if no out directory is specified the VGs are stored in
                # the current directory
                if args.out == "": 
                    args.out = os.path.abspath("./")

                workflow: BuildVG = BuildVG(args)

                if args.verbose:
                    end_args_parse: float = time.time()
                    print("Arguments parsed in %.2fs." % 
                          (end_args_parse - start_args_parse))
            # end if
        # end if
        
        #---------------------- findmotif options -----------------------#

        findmotif_err_msg: str = "Unexpected arguments for \"grafimo findmotif\": \"{}\""

        if args.workflow == "findmotif":
            if args.linear_genome:
                parser.error(findmotif_err_msg.format("-l, --linear-genome"))
                die(1)
            elif args.vcf:
                parser.error(findmotif_err_msg.format("-v, --vcf"))
                die(1)
            elif args.chroms_build:
                parser.error(findmotif_err_msg.format("--chroms-build"))
            elif args.chroms_prefix_build: 
                parser.error(findmotif_err_msg.format("--chroms-prefix-build"))
            elif args.chroms_namemap_build != NOMAP:
                parser.error(findmotif_err_msg.format("--chroms-namemap-build"))
            elif args.reindex:  # if default ignored
                parser.error(findmotif_err_msg.format("--reindex"))
                die(1)
            elif not args.graph_genome_dir and not args.graph_genome:
                parser.error(
                    "No arguments given for both \"--genome-graph\" and \"--genome-graph-dir\""
                )
                die(1)
            elif not args.bedfile:
                parser.error("No BED file given")
                die(1)
            elif not args.motif:
                parser.error("No motif PWM given")
                die(1)
            else:
                # only one between graph_genome and graph_genome_dir is allowed
                if args.graph_genome and args.graph_genome_dir:
                    parser.error(
                        "Only one argument between \"--genome-graph\" and \"--genome-graph-dir\""
                        " can be used"
                    )
                    die(1)

                # genome graph
                if args.graph_genome:
                    if (args.graph_genome.split('.')[-1] != "xg" and
                            args.graph_genome.split('.')[-1] != "vg"):
                        parser.error(
                            "Unrecognized genome variation graph format. Only"
                            "VG and XG format are allowed"
                        )
                        die(1)
                    elif not os.path.isfile(args.graph_genome):
                        parser.error(
                            "Unable to locate {}".format(args.graph_genome)
                        )
                        die(1)
                    else:
                        # using absolute path avoid potential problems
                        args.graph_genome = os.path.abspath(args.graph_genome)

                # genome graphs directory
                if args.graph_genome_dir:
                    if not os.path.isdir(args.graph_genome_dir):
                        parser.error(
                            "Unable to locate {}".format(args.graph_genome_dir)
                        )
                        die(1)
                    if len(glob(os.path.join(args.graph_genome_dir, "*.xg"))) <= 0:
                        parser.error(
                            "No genome variation graph found in {}".format(args.graph_genome_dir)
                        )
                        die(1)
                    else:
                        # using absolute path avoid potential problems
                        args.graph_genome_dir = os.path.abspath(args.graph_genome_dir)

                # BED file
                if args.bedfile:
                    if not isbed(args.bedfile, args.debug):
                        parser.error(
                            "The genomic coordinates must be given in UCSC BED files"
                        )
                        die(1)
                    else:
                        if not os.path.isfile(args.bedfile):
                            parser.error(
                                "Unable to locate {}".format(args.bedfile)
                            )
                else:
                    parser.error("No BED file given")

                # motif pwm
                if not args.motif:
                    parser.error("No motif PWM given")

                else:
                    motifs: List[str] = args.motif
                    for m in motifs:
                        if not isMEME_ff(m, args.debug) and not isJaspar_ff(m, args.debug):
                            parser.error(
                                "Unrecognized motif PWM file format. "
                                "{} does not follow the MEME or JASPAR format rules".format(m)
                            )
                            die(1)
                        if not os.path.isfile(m):
                            parser.error("Unable to locate {}".format(m))

                # background file
                if args.bgfile != UNIF:
                    if not os.path.isfile(args.bgfile):
                        parser.error("Unable to locate {}".format(args.bgfile))

                # pseudocount
                if args.pseudo <= 0:
                    parser.error(
                        "Pseudocount values must be > 0, got {}".format(args.pseudo))
                    die(1)

                # statistical significance threshold
                if args.threshold <= 0 or args.threshold > 1:
                    parser.error(
                        "Motif statistical significance threshold must be between 0 and 1")
                    die(1)

                # q-value flag
                if (not isinstance(args.no_qvalue, bool) or
                        (args.no_qvalue != False and args.no_qvalue != True)):
                    parser.error("\"--qvalue\" accepts only True or False values")
                    die(1)

                # no reverse flag
                if (not isinstance(args.no_reverse, bool) or
                        (args.no_reverse != False and args.no_reverse != True)):
                    parser.error("\"--no-reverse\" accepts only True or False values")
                    die(1)

                # text only flag
                if (not isinstance(args.text_only, bool) or
                        (args.text_only != False and args.text_only != True)):
                    parser.error("\"--text-only\" accepts only True or False values")
                    die(1)

                # chromosome to consider during VG scan
                if len(args.chroms_find) == 0:
                    args.chroms_find = [ALL_CHROMS]  # use all chromosome 
                else:
                    if anydup(args.chroms_find):
                        parser.error(
                            "Duplicated chromosome names given to \"--chroms-find\""
                        )

                # chromosome name-map
                if args.chroms_namemap_find != NOMAP:
                    if not os.path.isfile(args.chroms_namemap_find):
                        parser.error(
                            "Unable to locate {}".format(args.chroms_namemap_find)
                        )
                if (args.chroms_prefix_find and args.chroms_namemap_find != NOMAP):
                    parser.error(
                        "\"--chroms-prefix-find\" and \"chroms-namemap-find\" "
                        "cannot used together. Choose one of those options"
                    )

                # recomb flag
                if (not isinstance(args.recomb, bool) or
                        (args.recomb != False and args.recomb != True)):
                    parser.error("\"--recomb\" accepts only True or False values")
                    die(1)

                # out directory
                if args.out == "":  # default option
                    args.out = DEFAULT_OUTDIR 
                    print(args.out)

                # threshold on q-value flag
                if (not isinstance(args.qval_t, bool) or
                        (args.qval_t != False and args.qval_t != True)):
                    parser.error("\"--qvalueT accepts only True or False values")
                    die(1)
                elif args.no_qvalue == True and args.qval_t == True:
                    parser.error(
                        "Unable to apply statistical significance threshold on"
                        " q-values if you don't want them")
                    die(1)

                # number of graph regions to store as PNG images
                if args.top_graphs < 0:
                    parser.error("Negative number of regions to display")
                
                workflow: Findmotif = Findmotif(args)

                if args.verbose:
                    end_args_parse: float = time.time()
                    print("Arguments parsed in %.2fs." % 
                          (end_args_parse - start_args_parse))
            # end if
        # end if

        # chck that external dependencies are satisfied
        if args.verbose:
            sys.stderr.write("Checking GRAFIMO external dependencies {}\n".format(EXT_DEPS))
            start_deps: float = time.time()
        satisfied: bool 
        deps_lack: List[str] 
        satisfied, deps_lack = check_deps()
        if not satisfied and len(deps_lack) > 0:
            errmsg = "Some dependencies are not satisfied: {}.\nPlease solve them before running GRAFIMO.\n"
            exception_handler(DependencyError, errmsg.format(deps_lack), args.debug)
        elif not satisfied and len(deps_lack) <= 0:
            errmsg = "Dependencies satisfied, but unable to recover them.\n Be sure they are in system PATH.\n"
            exception_handler(DependencyError, errmsg, args.debug)

        if args.verbose and satisfied:
            end_deps: float = time.time()
            print("Dependencies satisfied.")
            print("Dependencies checked in %.2fs." % (end_deps - start_deps))

        #---------------------------------------------------------------
        # dependency check was ok, so we go to workflow selection:
        #   * construction of the genome variation graph for 
        #     each chromosome or a user defined subset of them
        #   * scan of a precomputed VG or a set of precomputed VG
        if isinstance(workflow, BuildVG): buildvg(workflow, args.debug)
        elif isinstance(workflow, Findmotif): findmotif(workflow, args.debug)
        else:
            errmsg = "Expected BuildVG or Findmotif, got {}.\n"
            exception_handler(TypeError, errmsg.format(type(workflow).__name__), args.debug)
        
        end: float = time.time()  # GRAFIMO execution finishes here
        print("Elapsed time %.2fs." % (end - start))

    except KeyboardInterrupt:
        sigint_handler()
    finally:
        pass
# end of main()


#-----------------------------------------------------------------------
# Entry point for GRAFIMO
#

if __name__ == "__main__":
    main()

