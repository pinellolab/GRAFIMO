"""

GRAFIMO version {version}

Copyright (C) 2019 Manuel Tognon <manu.tognon@gmail.com> <manuel.tognon@studenti.univr.it>

GRAph-based Find Individual Motif Occurrences scores sequences likely to be 
transcription factor binding sites in ChIP-seq narrow peaks

Usage:
    grafimo --linear-genome LINEAR-GENOME --vcf VCF --bedfile BEDFILE --motif MOTIF [options]
    
    grafimo --graph-genome GRAPH-GENOME --bedfile BEDFILE --motif MOTIF [options]

    grafimo --graph-genome-dir PATH-TO-GRAPH-GENOME --bedfile BEDFILE --motif MOTIF [options]
    
The tool takes in input a linear genome (in FA or FASTA format) and a VCF file
(must be compressed, e.g. vcf-file.vcf.gz) to construct your own graph genomes,
or a complete graph_genome (better if in XG format), or a set of graph genomes (one for each
chromosome); the latter option is suggested to run GRAFIMO on your laptop.
The results will be written in the directory specified by the user
(--out option or default value 'grafimo_out') and they are made of a TSV file, a
GFF file and an HTML version of the TSV.

The tool is also error tolerant in the subgraphs extraction, since it doesn't
stop its execution although we can have some exceptions.

Citation:
    
    xxxx. Xxxx, xxxx
    
    
Run 'grafimo --help'to see all command-line options.
See https://github.com/pinellolab/GRAFIMO or https://github.com/InfOmics/GRAFIMO for the full documentation.

"""

from argparse import ArgumentParser, SUPPRESS, HelpFormatter
from grafimo.grafimo import __version__, with_vg_pipeline, without_vg_pipeline
from grafimo.GRAFIMOException import DependencyException, VGException, PipelineException
from grafimo.utils import die, check_deps, EXT_DEPS
import time
import sys
import os
import glob
import multiprocessing as mp
import grafimo.vgCreation as vgc


class GRAFIMOArgumentParser(ArgumentParser):
    """
    
        This redefinition of the arguments parser:
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

def get_AP():
    """
        Return the parser containing the input arguments
        ----
        Parameters:
            None
        ----
        Returns:
            parser (GRAFIMOArgumentParser)
    
    """
    
    parser = GRAFIMOArgumentParser(usage=__doc__, add_help=False)
   
    group = parser.add_argument_group("Options")
    group.add_argument("-h", "--help", action="help", help="Show the help message and exit")
    group.add_argument("--version", action="version", help="Show the version and exit", version=__version__)
    group.add_argument("--cores", type=int, default=0, metavar="NCORES", nargs='?', const=0,
                           help="Number of cores to use. The default value is %(default)s (auto-detection)."
                                "If you chose to query the whole genome variation graph note that the default "
                                "option is to use only one core to avoid memory issues; if you are running "
                                "GRAFIMO on a machine with RAM >= 128 GB maybe you can use more cores")
    
    # the following arguments depend on the user approach
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
    
    # optional arguments
    group.add_argument("-k", "--bgfile", type=str,
                        help='Path to a file which defines the background distribution [optional]', nargs='?',
                        const='', default='', metavar='BACKGROUND')
    group.add_argument("-p", "--pseudo", type=float,
                       help='Pseudocount to add to the motif counts [optional]', nargs='?', const='0.1',
                       default='0.1', metavar='PSEUDOCOUNT')
    group.add_argument("-t", "--pvalueT", type=float, nargs='?', default=1e-4,
                        metavar='THRESHOLD', const='1e-4',
                        help='Defines a threshold on the p-value obtained form the scoring ' 
                                'of each sequence.\n'
                                'Default is 1e-4 [optional]')
    group.add_argument("-q", "--no-qvalue", action='store_true', default=False, dest='no_qvalue',
                        help='Flag parameter, if it is set to True will be computed the q-value for each hit')
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

    group.add_argument("--top-graphs", type=int,
                        help='Number of graphs of the regions that will be stored as PNG images in the results '
                        'directory [optional]',
                        nargs='?', const=0, default=0, metavar='GRAPHS_NUM', dest='top_graphs')

    group.add_argument("--verbose",  default = False, action = 'store_true',
                        help = "Output a lot of additional informations about the execution")
    
    return parser

def main(cmdLineargs = None):
    """
    
        Main function of the tool
        ----
        Parameters:
            cmdLineArgs (str) : the arguments given in command line 
        ----
        Returns:
            None
    
    """

    start = time.time() # compute elapsed time

    WITH_VG_CREATION = False
    
    parser = get_AP()
    if cmdLineargs is None:
        cmdLineargs = sys.argv[1:] # take input args
        
    args = parser.parse_args(args = cmdLineargs)
    
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
        cores = mp.cpu_count() # take all the CPUs available

    elif args.cores == 0 and args.graph_genome:
        cores = 1 # the default if we make a query on the whole genome is 1 core
                  # this slow doen a lot the query process but is useful to avoid
                  # memory problems
        
    else:
        cores=args.cores
        
    if args.linear_genome:
        if  (args.linear_genome.split('.')[-1] != 'fa' and \
            args.linear_genome.split('.')[-1] != 'fasta'):
            parser.error('The linear genome must be in fasta format (.fasta and .fa allowed)')
            
        else:
            linear_genome = args.linear_genome

            if len(glob.glob(linear_genome)) != 1:
                parser.error('Cannot find the specified linear genome file')
            
    if args.chroms:
        chroms = args.chroms
        
    else:
        chroms = []
            
    if args.vcf:
        if args.vcf.split('.')[-1] != 'gz' or \
            args.vcf.split('.')[-2] != 'vcf': # allow only compressed VCF files
            parser.error('Incorrect vcf file given')
        else:
            vcf = args.vcf

            if len(glob.glob(vcf)) <= 0:
                parser.error('Cannot find the specified VCF file')
            
            WITH_VG_CREATION = True
            
    if args.graph_genome:
        if (args.graph_genome.split('.')[-1] != 'xg' and \
            args.graph_genome.split('.')[-1] != 'vg'):
            parser.error("Incorrect genome graph specified. Please check that they are in XG or VG format")

        elif (not os.path.isfile(args.graph_genome)):
            parser.error("Unable to find the specified graph genome")
        
        else:
            graph_genome = args.graph_genome
            WITH_VG_CREATION = False
            
    if args.graph_genome_dir:
        if (not os.path.isdir(args.graph_genome_dir)):
            parser.error("The specified directory where to look for graph genomes doesn't exist")

        if args.graph_genome_dir[-1] == '/':
            graph_genome_dir = args.graph_genome_dir

        else:
            graph_genome_dir = ''.join([args.graph_genome_dir, '/'])
            
        if len(glob.glob(graph_genome_dir+'*.xg')) <= 0:
            parser.error('No xg file found in the specified directory')
        
        WITH_VG_CREATION = False # in any case we skip the VG creation step
        
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
                parser.error('Cannot motif file: ' + m)
        
    if args.bgfile:
        bgfile = args.bgfile # we have a path to the bg file

        if len(glob.glob(bgfile)) <= 0:
            parser.error('Cannot find the specified background file')

    else:
        bgfile = args.bgfile # empty string
        
    if args.pseudo <= 0:
        parser.error('The pseudocount cannot be less than or equal 0')
    
    else:
        pseudocount = args.pseudo
        
    if args.pvalueT <= 0 or args.pvalueT > 1:
        parser.error('The pvalue threshold must be between 0 and 1')
    
    else:
        pvalueT = args.pvalueT
        
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
        # to the name
        dest = ''.join([dest, '_', str(os.getpid())])

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
                                    "\nPlease solve them before running GRAFIMO")
        die(1)

    elif not sat and len(deps_lack) <= 0:
        raise DependencyException("Some dependencies were not found, but was not possible to track them")
        die(1)

    if verbose and sat:
        print("Dependencies correctly satisfied")

    # the dependency check was OK

    if WITH_VG_CREATION:
        
        if verbose:
            print("\nEntering the pipeline with the variation graph creation\n")

        with_vg_pipeline(cores, linear_genome, vcf, chroms, bedfile, motifs, bgfile,
                            pseudocount, pvalueT, no_reverse, qvalue, text_only, dest,
                            top_graphs, WITH_VG_CREATION, verbose)
        
        
    elif not WITH_VG_CREATION:

        if verbose:
            print("\nEntering the pipeline without the variation graph creation\n")
        
        if args.graph_genome:
        
            if graph_genome.split('.')[-1] != 'xg' and \
                graph_genome.split('.')[-1] != 'vg':
                raise VGException("The genome graph must be in VG or XG format")
                die(1)
        
            elif graph_genome.split('.')[-1] == 'vg': # we are given a vg genome
                vg = graph_genome
                xg = vgc.indexVG(vg)
            
            else: # we are given an xg genome
                xg = graph_genome

            if verbose:
                print("The graph " + xg + " will be queried\n")
        
            without_vg_pipeline(cores, xg, bedfile, motifs, bgfile, pseudocount, 
                                    pvalueT, no_reverse, qvalue, text_only, dest, top_graphs, WITH_VG_CREATION)
            
        elif args.graph_genome_dir:
            
            gplus=True # defines if the input is a single xg or more than one

            if verbose:
                print("The graphs contained in directory " + graph_genome_dir + " will be queried\n")
            
            without_vg_pipeline(cores, graph_genome_dir, bedfile, motifs, bgfile, pseudocount,
                                    pvalueT, no_reverse, qvalue, text_only, dest, top_graphs, WITH_VG_CREATION, gplus, chroms)
            
    else:
        # error in the pipeline flag
        msg="\n\nWrong input. Unable to determine which pipeline to follow"
        raise PipelineException(msg)
        die(1)
        
    end=time.time()
    
    print('\nelapsed time', end-start)
    
    
###########################################
#
# Entry point for GRAFIMO
#
###########################################

if __name__=="__main__":
    main()
    
