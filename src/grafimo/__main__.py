#!/usr/bin/env python
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

GRAph-based Find Individual Motif Occurrences scores sequences likely to be 
transcription factor binding sites from ChIP-seq narrow peaks data.

Usage:
    grafimo --linear-genome LINEAR-GENOME --vcf VCF --bedfile BEDFILE --motif MOTIF [options]

    grafimo --graph-genome-dir PATH-TO-GRAPH-GENOME --bedfile BEDFILE --motif MOTIF [options]

    grafimo --graph-genome GRAPH-GENOME --bedfile BEDFILE --motif MOTIF [options]
    
Construct the variation graph from a linear reference genome and your own VCF file.
The genome is splitted in the 23 chromosomes (a graph for each chromosome will be built).
From the genome graphs get all the sequences contained in the regions defined in your 
BED file.
If you have a pre-built genome graph, the scanning can be done on it; the scanning can also 
be done on a set of genome graphs, assuming that each genome graph represent a single chromosome,
by giving the directory were they are stored. Note that the latter approach (that is the one used
by the full pipeline) has shown the best performances.
The results will be written in a directory, named grafimo_out_PID_TFCode; in the directory you'll
find three files: 
    - grafimo_out.csv
    - grafimo_out.html
    - grafimo_out.gff

Citation:
    
    xxxx. Xxxx, xxxx
    
    
Run 'grafimo --help'to see all command-line options.

See https://github.com/pinellolab/GRAFIMO or https://github.com/InfOmics/GRAFIMO for the full documentation.

"""

from argparse import ArgumentParser, SUPPRESS, HelpFormatter
from grafimo.grafimo import __version__, with_vg_pipeline, without_vg_pipeline
from grafimo.GRAFIMOException import DependencyException, VGException, PipelineException
from grafimo.utils import die, check_deps, initialize_chroms_list, EXT_DEPS
import time
import sys
import os
import glob
import multiprocessing as mp
import grafimo.vgCreation as vgc


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
                
    def __init__(self, *args, **kwargs):
        kwargs['formatter_class'] = self.GRAFIMOHelpFormatter
        kwargs['usage'] = kwargs['usage'].replace("{version}", __version__)
        super().__init__(*args, **kwargs)
        
    def error(self, msg):
        sys.stderr.write("Run 'grafimo --help' to see the usage\n")
        self.exit(2, "\n{0}: ERROR: {1}\n".format(self.prog, msg))
        

class CmdError(Exception):
    pass


class Pipeline(object):
    """
        Generic pipeline object
    """
    def __init__(self):
        pass


class withVGCreationPipeline(Pipeline):
    """
        Object to represent the pipeline with 
        the creation of the genome graphs
    """
    pass


class onlyMotifScanningPipeline(Pipeline):
    """
        Object to represent the pipeline that
        does only the motif scanning on the given 
        genome graph(s)
    """
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
    
    # define the parser
    parser = GRAFIMOArgumentParser(usage=__doc__, add_help=False)
   
    #####################################################################
    # start reading the arguments
    #####################################################################
    group = parser.add_argument_group("Options")
    group.add_argument("-h", "--help", action="help", help="Show the help message and exit")
    group.add_argument("--version", action="version", help="Show the version and exit", version=__version__)
    group.add_argument("--cores", type=int, default=0, metavar="NCORES", nargs='?', const=0,
                           help="Number of cores to use. The default value is %(default)s (auto-detection)."
                                "If you chose to query the whole genome variation graph note that the default "
                                "option is to use only one core to avoid memory issues; if you are running "
                                "GRAFIMO on a machine with RAM >= 128 GB maybe you can use more cores")
    
    #######################
    # mandatory arguments
    #######################

    # the following arguments depend on the user approach (only motif scanning || genome graph creation + motif scanning)
    group.add_argument("-l", "--linear-genome", type=str, help='Path to the linear genome (fasta format required)',
                           nargs='?', default='', metavar='LINEAR-GENOME', dest='linear_genome')
    group.add_argument("-c", "--chroms", type=str, nargs='*',
                        help='Chromosome numbers to consider in the pipeline application.\n'
                             'To consider all the chromosomes, just skip this argument',
                        default=[], metavar=('1', 'X'))
    group.add_argument("-v", "--vcf", type=str, help='Path to the VCF file. The vcf must be compressed (vcf.GZ)',
                        nargs='?', default='', metavar='VCF')
    group.add_argument("-g", "--graph-genome", type=str, nargs= '?', metavar='GRAPH-GENOME',
                           help='Path to the VG genome graph file (vg or xg format required)',
                           default='', dest='graph_genome')
    group.add_argument("-d", "--graph-genome-dir", type=str, nargs='?', metavar='GRAPH-GENOMES-DIR',
                           help='Path to a directory containing a variable number of VG ' 
                                'graph genomes (only XG format is allowed)',
                           default='', dest='graph_genome_dir')
    group.add_argument("-b", "--bedfile", type=str, help='Path to the BED file that defines the regions to score',
                           metavar='BEDFILE')
    group.add_argument("-m", "--motif", type=str, nargs='+', metavar=('MOTIF1', 'MOTIF2'),
                           help='Path to the motif file to use for the scoring of the '
                                'sequences (JASPAR or MEME format required)')
    
    ######################
    # optional arguments
    ######################
    group.add_argument("-k", "--bgfile", type=str,
                        help='Path to a file which defines the background distribution [optional]', nargs='?',
                        const='', default='', metavar='BACKGROUND')
    group.add_argument("-p", "--pseudo", type=float,
                       help='Pseudocount to add to the motif counts [optional]', nargs='?', const='0.1',
                       default='0.1', metavar='PSEUDOCOUNT')
    group.add_argument("-t", "--threshold", type=float, nargs='?', default=1e-4,
                        metavar='THRESHOLD', const='1e-4',
                        help='Defines a threshold on the p-value (by default) obtained form the scoring ' 
                                'of each sequence. It is possible to apply the threshold on the q-value '
                                'using the -q (--qvalueT) option.\n'
                                'Default is 1e-4 [optional]')
    group.add_argument("-q", "--no-qvalue", action='store_true', default=False, dest='no_qvalue',
                        help='Flag parameter, if it is set to True the q-value for each hit will not be computed')
    group.add_argument("-r", "--no-reverse", default=False, action='store_true', dest='no_reverse',
                        help='Flag parameter, if it is set to True will be scored only '
                                'sequences from the forward strand')
    group.add_argument("-f", "--text-only", default=False, action='store_true',
                       dest='text_only',
                       help="Print the results only on the terminal "
                            "without writing the results files",
                       )
    group.add_argument("-o", "--out", type=str,
                        help='Name of the directory where the results will be stored [optional]',
                        nargs='?', const='grafimo_out', default='grafimo_out', metavar='OUTDIR')
    group.add_argument("--qvalueT", action='store_true', default=False, dest='qval_t',
                       help='Flag parameter, if set to true the threshold will be applied on the '
                            'q-value, rather than the p-value ')
    group.add_argument("--top-graphs", type=int,
                        help='Number of graphs of the regions that will be stored as PNG images in the results '
                        'directory [optional]',
                        nargs='?', const=0, default=0, metavar='GRAPHS_NUM', dest='top_graphs')

    group.add_argument("--verbose",  default = False, action = 'store_true',
                        help = "Output a lot of additional informations about the execution")
    
    return parser


def guess_pipeline_from_args(args):
    """
        Return the pipeline to follow from the input
        arguments
        ----
        Parameters:
            - args (GRAFIMOArgumentParser)
        ----
        Return:
            - pipeline (str)
    """

    pipeline = None

    if args.linear_genome and args.vcf:
        pipeline = withVGCreationPipeline()
    
    if args.graph_genome_dir or args.graph_genome:
        pipeline = onlyMotifScanningPipeline()

    return pipeline


def main(cmdLineargs = None):
    """
    
        Main function of GRAFIMO.

        The arguments given in input are checked for consistency,
        then a pipeline is followed.

        ----
        Parameters:
            cmdLineArgs (str) 
        ----
        Returns:
            None
    
    """

    start = time.time() # the run begin here
    
    # read command-line arguments
    parser = get_AP()

    if cmdLineargs is None:
        cmdLineargs = sys.argv[1:] # take input args
        
    args = parser.parse_args(args = cmdLineargs)
    
    #####################################################
    # checking arguments consistency
    #####################################################

    # and a long series of if ..  else begins
    # good luck !!!
    if not args.linear_genome and not args.graph_genome and not args.graph_genome_dir:
        parser.error("Needed at least one between '--linear_genome --vcf [other args]', '--graph_genome [other args]' and '--graph_genome_dir [other args]'")
        
    if args.linear_genome and args.graph_genome:
        parser.error("Only one between --linear_genome --vcf [other args]' and '--graph_genome [other args]' is allowed")
        
    if args.linear_genome and args.graph_genome_dir:
        parser.error("Only one between --linear_genome --vcf [other args]' and '--graph_genome_dir [other args]' is allowed")
        
    if args.graph_genome and args.vcf:
        parser.error("Only one between --linear_genome --vcf [other args]' and '--graph_genome [other args]' is allowed")
    
    if args.graph_genome_dir and args.vcf:
        parser.error("Only one between --linear_genome --vcf [other args]' and '--graph_genome_dir [other args]' is allowed")
        
    if args.graph_genome and args.graph_genome_dir:
        parser.error("Only one between '--graph_genome [other args]' and '--graph_genome_dir [other args]'is allowed")
        
    if args.chroms and args.graph_genome:
        parser.error("The --chroms option can be used only when creating the genome graph or querying the graph of single chromosomes")

    if args.linear_genome and not args.vcf:
        parser.error('No VCF file given')
        
    if args.cores < 0:
        parser.error('The number of cores cannot be negative')
        
    if args.cores == 0:
        cores = mp.cpu_count() # take all the available CPUs

    elif args.cores == 0 and args.graph_genome:
        cores = 1 # the default, if we make work with a whole genome-graph is 1 core
                  # this slows down a lot the region query process but we are sure to avoid
                  # memory problems
                  # 
                  # ****
                  # If you call more cores querying the whole genome graph 
                  # BE SURE TO HAVE > 128 GB of RAM !!!!!!
                  # ****
        
    else:
        cores = args.cores
        
    if args.linear_genome:
        if  (args.linear_genome.split('.')[-1] != 'fa' and \
            args.linear_genome.split('.')[-1] != 'fasta'):
            parser.error('The linear genome must be in fasta format (.fasta and .fa allowed)')
            
        else:
            linear_genome = args.linear_genome

            if len(glob.glob(linear_genome)) != 1:
                parser.error('Cannot find the specified linear genome file')

            # it's safer to get the absolute path to the linear genome
            linear_genome = os.path.abspath(linear_genome)
            
    # working with the chromosomes list
    chroms = initialize_chroms_list(args.chroms)
            
    if args.vcf:
        if args.vcf.split('.')[-1] != 'gz' or \
            args.vcf.split('.')[-2] != 'vcf': # allow only compressed VCF files
            parser.error('Incorrect vcf file given')
        else:
            vcf = args.vcf

            if len(glob.glob(vcf)) <= 0:
                parser.error('Cannot find the specified VCF file')

            # it's safer to get the absolute path to the VCF
            vcf = os.path.abspath(vcf)
            
    if args.graph_genome:
        if (args.graph_genome.split('.')[-1] != 'xg' and \
            args.graph_genome.split('.')[-1] != 'vg'):
            parser.error("Incorrect genome graph specified. Please check that they are in XG or VG format")

        elif (not os.path.isfile(args.graph_genome)):
            parser.error("Unable to find the specified graph genome")
        
        else:
            graph_genome = args.graph_genome
            graph_genome = os.path.abspath(graph_genome)
            
    if args.graph_genome_dir:
        if (not os.path.isdir(args.graph_genome_dir)):
            parser.error("The specified directory where to look for graph genomes doesn't exist")

        if args.graph_genome_dir[-1] == '/':
            graph_genome_dir = args.graph_genome_dir

        else:
            graph_genome_dir = ''.join([args.graph_genome_dir, '/'])
            
        if len(glob.glob(graph_genome_dir+'*.xg')) <= 0:
            parser.error('No xg file found in the specified directory')

        graph_genome_dir = os.path.abspath(graph_genome_dir)
        
    if args.bedfile:
        if args.bedfile.split('.')[-1] != 'bed':
            parser.error('Incorrect bedfile given')
    
        else:
            bedfile = args.bedfile

            if len(glob.glob(bedfile)) <= 0:
                parser.error('Cannot find the specified BED file')
    
    else:
        parser.error('No BED file given')
        
    if not args.motif:
        parser.error('No motif given')
        
    else:
        motifs = args.motif

        # check if the given motifs exist
        for m in motifs:
            if len(glob.glob(m)) <= 0:
                parser.error('Cannot find motif file: ' + m)
        
    if args.bgfile:
        bgfile = args.bgfile # we have a path to the bg file

        if len(glob.glob(bgfile)) <= 0:
            parser.error('Cannot find the specified background file')

    else:
        bgfile = args.bgfile # empty string => uniform distribution
        
    if args.pseudo <= 0:
        parser.error('The pseudocount cannot be less than or equal 0')
    
    else:
        pseudocount = args.pseudo
        
    if args.threshold <= 0 or args.threshold > 1:
        parser.error('The pvalue threshold must be between 0 and 1')
    
    else:
        threshold = args.threshold
        
    if not isinstance(args.no_qvalue, bool) and args.no_qvalue != False \
        and args.no_qvalue != True:
        parser.error('The --qvalue parameter accepts only True or False as values')
        
    else:
        qvalue=not(bool(args.no_qvalue))
        
    if not isinstance(args.no_reverse, bool) or (args.no_reverse != False \
        and args.no_reverse != True):
        parser.error('The --no-reverse parameter accepts only True or False as values')
        
    else:
        no_reverse = bool(args.no_reverse)

    if not isinstance(args.text_only, bool) and args.text_only != False \
        and args.text_only != True:
        parser.error('The --text-only parameter accepts only True or False values')

    else:
        text_only = bool(args.text_only)
        
    if args.out:
        dest = args.out

    if dest == 'grafimo_out': # default option

        # to make unique the output directory we add the PID
        # to the name.
        #
        # This is useful when calling grafimo in different runs on the 
        # same machine.
        #
        # believe me, this can happen

        dest = ''.join([dest, '_', str(os.getpid())])

    if not isinstance(args.qval_t, bool) or(args.qval_t != False \
        and args.qval_t != True):
        parser.error("The --qvalueT parameter accepts only True or False as values")

    elif args.no_qvalue == True and args.qval_t == True:
        parser.error("Cannot apply the threshold on the q-values if you don't want them")

    else:
        qval_t = bool(args.qval_t)

    if args.top_graphs < 0:
        parser.error("The number of region graphs to show must be positive")

    else:
        top_graphs = args.top_graphs

    if not isinstance(args.verbose, bool) and args.verbose != False \
        and args.verbose != True:
        parser.error('The --verbose parameter accepts only True or False values')

    else:
        verbose = bool(args.verbose)

    # checking that external dependencies are satisfied 
    if verbose:
        print("Checking GRAFIMO external dependencies " + str(EXT_DEPS))

    sat, deps_lack = check_deps()
    
    if not sat and len(deps_lack) > 0:
        raise DependencyException("\n\nERROR: The following dependencies are not sastisfied: " + 
                                    str(deps_lack) + 
                                    "\nPlease, solve them before running GRAFIMO")
        die(1)

    elif not sat and len(deps_lack) <= 0:
        raise DependencyException("Some dependencies were not found, but was not possible to track them")
        die(1)

    if verbose and sat:
        print("Dependencies correctly satisfied")

    # the dependency check was OK

    # choose the pipeline from the arguments
    pipeline = guess_pipeline_from_args(args)

    if pipeline == None:
        raise PipelineException("Don't know which pipeline to follow")

    if isinstance(pipeline, withVGCreationPipeline):
        if verbose:
            print("\nEntering the pipeline with the variation graph creation\n")

        with_vg_pipeline(cores, linear_genome, vcf, chroms, bedfile, motifs, bgfile,
                            pseudocount, threshold, no_reverse, qvalue, qval_t, text_only, dest,
                            top_graphs, True, verbose)
        
    elif isinstance(pipeline, onlyMotifScanningPipeline):

        if verbose:
            print("\nEntering the pipeline without the variation graph creation\n")
        
        if args.graph_genome:
        
            if graph_genome.split('.')[-1] != 'xg' and \
                graph_genome.split('.')[-1] != 'vg':
                raise VGException("The genome graph must be in VG or XG format")
                die(1)
        
            elif graph_genome.split('.')[-1] == 'vg': 
                # we are given a vg genome, then we have to index it
                vg = graph_genome
                vgc.indexVG(vg)
            
            else: # we are given an xg genome
                xg = graph_genome

            if verbose:
                print("The graph " + xg + " will be queried\n")
        
            without_vg_pipeline(cores, xg, bedfile, motifs, bgfile, pseudocount, 
                                    threshold, no_reverse, qvalue, qval_t, text_only, dest, top_graphs, False)
            
        elif args.graph_genome_dir:
            
            gplus = True # defines if the input is a single xg or more than one

            if verbose:
                print("The graphs contained in directory " + graph_genome_dir + " will be queried\n")
            
            without_vg_pipeline(cores, graph_genome_dir, bedfile, motifs, bgfile, pseudocount,
                                    threshold, no_reverse, qvalue, qval_t, text_only, dest, top_graphs, False,
                                    gplus, chroms)
            
    else:
        # error in the pipeline flag
        msg = "\n\nWrong input. Unable to determine which pipeline to follow"
        raise PipelineException(msg)
        die(1)
        
    end = time.time() # the run  ends here
    
    print('\nelapsed time', end-start)
    
    
###########################################
#
# Entry point for GRAFIMO
#
###########################################

if __name__=="__main__":
    main()
    
