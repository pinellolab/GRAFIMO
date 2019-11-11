"""

GRAFIMO version {version}

Copyright (C) 2019 Manuel Tognon <manu.tognon@gmail.com> <manuel.tognon@studenti.univr.it>

GRAph-based Find Individual Motif Occurrences scores sequences likely to be 
transcription binding sites in ChIP-seq peaks

Usage:
    grafimo --cores [NCORES] --linear_genome [LINEAR_GENOME.fa] --vcf [VCF.vcf.gz] --bedfile [BEDFILE.bed] --motif [MOTIF.jaspar] [options]
    
    grafimo --cores [NCORES] --graph_genome [GRAPH_GENOME.xg] --bedfile [BEDFILE.bed] --motif [MOTIF.jaspar] [options]
    
The tool takes in input both a linear genome (in .fa format) and a vcf file 
(must be compressed, e.g. vcffile.vcf.gz) or a graph_genome (better if in .xg).
The results will be written to a .csv file in the directory specified by the user
(--o option or defalut value 'grafimo_out').

The tool is also error tollerant in the subgrapahs extraction, since it doesn't 
stop its execution although we can have some exceptions.

Citation:
    
    xxxx. Xxxx, xxxx
    
    
Run 'grafimo --help'to see all command-line options.
See https://github.com/InfOmics/GRAFIMO for the full documentation.

"""

from argparse import ArgumentParser, SUPPRESS, HelpFormatter
from grafimo.grafimo import __version__, with_vg_pipeline, without_vg_pipeline
from grafimo.GRAFIMOException import VGException, PipelineException
from grafimo.utils import die
import time
import sys
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
                args=usage, actions, groups, ''
                self._add_item(self._format_usage, args)
                
    def __init__(self, *args, **kwargs):
        kwargs['formatter_class']=self.GRAFIMOHelpFormatter
        kwargs['usage']=kwargs['usage'].replace("{version}", __version__)
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
            parser(GRAFIMOArgumentParser)
    
    """
    
    parser=GRAFIMOArgumentParser(usage=__doc__, add_help=False)
   
    group=parser.add_argument_group("Options")
    group.add_argument("-h", "--help", action="help", help="Show the help message and exit")
    group.add_argument("--version", action="version", help="Show the version and exit", version=__version__)
    group.add_argument("--cores", type=int, default=0, metavar='NCORES', nargs='?', const=0, 
                           help="Number of cores to use. The defalut value is 0 (auto-detection)")
    
    # the three following arguments depend on the user approach
    group.add_argument("-l", "--linear_genome", type=str, help='path to the linear genome (fasta format required)',
                           nargs='?', default='', metavar='LINEAR-GENOME.fa')
    group.add_argument("-c", "--chroms", type=str, help='chromosomes to consider in the pipeline application (to take all chromosomes argument not needed)',
                           nargs='*', default=[], metavar='1, 2, X, Y')
    group.add_argument("-v", "--vcf", type=str, help='path to your .vcf.gz file', nargs='?', default='', metavar='VCF.vcf.gz')
    group.add_argument("-g", "--graph_genome", type=str, nargs= '?', metavar='PATH-TO-GRAPH-GENOME.xg',
                           help='path to the VG graph genome file (vg or xg format required)',
                           default='') 
    group.add_argument("-d", "--graph_genome_dir", type=str, nargs='?', metavar='PATH-TO-GRAPH-GENOMES-DIRECTORY',
                           help='path to a directory containing a variable number of VG graph genomes (only XG format is allowed)',
                           default='')
    group.add_argument("-b", "--bedfile", type=str, help='path to the BED file that defines the regions to score',
                           metavar='BEDFILE.bed')
    group.add_argument("-m", "--motif", type=str, metavar='MOTIF_FILE',
                           help='path to the motif file to use for the scoring of the sequences (JASPAR format required)')
    
    # optional arguments
    group.add_argument("-k", "--bgfile", type=str, help='path to a file which defines the background distribution [optional]', nargs='?',
                          const='', default='', metavar='BGFILE')
    group.add_argument("-p", "--pseudo", type=float, help='pseudocount to add to the motif counts [optional]', nargs='?', const='0.1',
                       default='0.1', metavar='PSEUDOCOUNT')
    group.add_argument("-t", "--pvalueT", type=float, nargs='?', const='1e-3', metavar='PVALUE-THRESHOLD',
                           help='defines a threshold on the p-value obtained form the scoring of each sequence (1e-3 by default) [optional]',
                           default=1e-3)
    group.add_argument("-q", "--qvalue", type=bool, nargs='?', const=False, default=False, metavar='True | False',
                           help='flag parameter, if it is set to True will be computed the q-value for each hit')
    group.add_argument("-r", "--no_reverse", type=bool, nargs='?', const=False, default=False, metavar='True | False',
                           help='flag parameter, if it is set to True will be scored only sequences from the forward strand')
    group.add_argument("-o", "--out", type=str, help='name of the directory where the results will be stored [optional]',
                           nargs='?', const='grafimo_out', default='grafimo_out', metavar='OUTDIR')
    
    return parser

def main(cmdLineargs=None):
    """
    
        Main function of the tool
        ----
        Parameters:
            cmdLineArgs (str) : the arguments given in command line
            
        ----
        Returns:
            None
    
    """
    
    start=time.time()
    
    WITH_VG_CREATION=False
    
    parser=get_AP()
    if cmdLineargs is None:
        cmdLineargs=sys.argv[1:] # take input args
        
    args=parser.parse_args(args=cmdLineargs)
    
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
        parser.error("The --chroms option can be used only when creating the genome graph")
        
    if args.cores<0:
        parser.error('The number of cores cannot be negative')
        
    if args.cores==0:
        cores=mp.cpu_count() # take all the CPUs available
        
    else:
        cores=args.cores
        
    if args.linear_genome:
        if args.linear_genome.split('.')[-1]!='fa' and \
            args.linear_genome.split('.')[-1]!='fasta':
            parser.error('The linear genome must be in fasta format (.fasta and .fa allowed)')
            
        else:
            linear_genome=args.linear_genome
            
    if args.chroms:
        chroms=args.chroms
        
    else:
        chroms=[]
            
    if args.vcf:
        if args.vcf.split('.')[-1]!='gz' or \
            args.vcf.split('.')[-2]!='vcf': # allow also the compressed vcf files
                                            # like something.vcf.gz
            parser.error('Incorrect vcf file given')
        else:
            vcf=args.vcf
            WITH_VG_CREATION=True
            
            
    if args.graph_genome:
        if args.graph_genome.split('.')[-1]!='xg' and \
            args.graph_genome.split('.')[-1]!='vg':
            parser.error("Incorrect genome graph specified. Please check that they are in XG or VG format")
        
        else:
            graph_genome=args.graph_genome
            WITH_VG_CREATION=False
            
    if args.graph_genome_dir:
        if args.graph_genome_dir[-1]=='/':
            graph_genome_dir=args.graph_genome_dir

        else:
            graph_genome_dir=''.join([args.graph_genome_dir, '/'])
            
        WITH_VG_CREATION=False # in any case we skip the VG creation step
        
    if args.bedfile.split('.')[-1]!='bed':
        parser.error('Incorrect bedfile given')
    
    else:
        bedfile=args.bedfile
        
    if args.motif.split('.')[-1]!='jaspar' and \
        args.motif.split('.')[-1]!='meme':
        parser.error('Incorrect motif file given, only JASPAR or MEME format allowed')
        
    else:
        motif=args.motif
        
    if args.bgfile:
        bgfile=args.bgfile # we have a path to the bg file
    else:
        bgfile=args.bgfile # empty string
        
    if args.pseudo <= 0:
        parser.error('The pseudocount cannot be less than or equal 0')
    
    else:
        pseudocount=args.pseudo
        
    if args.pvalueT < 0 or args.pvalueT > 1:
        parser.error('The pvalue threshold must be between 0 and 1')
    
    else:
        pvalueT=args.pvalueT
        
    if not isinstance(args.qvalue, bool) and args.qvalue != 'False' \
        and args.qvalue != 'True':
        parser.error('The qvalue parameter accepts only True or False as values')
        
    else:
        qvalue=bool(args.qvalue)
        
    if not isinstance(args.no_reverse, bool) and args.no_reverse != 'False' \
        and args.no_reverse != 'True':
        parser.error('The no_reverse parameter accepts only True or False as values')
        
    else:
        no_reverse=bool(args.no_reverse)
        
    if args.out:
        dest=args.out
        
    if WITH_VG_CREATION:
            
        dest+='/' 
        
        with_vg_pipeline(cores, linear_genome, vcf, chroms, bedfile, motif, bgfile,
                             pseudocount, pvalueT, no_reverse, dest, WITH_VG_CREATION)
        
        
    elif not WITH_VG_CREATION:
        
        if args.graph_genome:
        
            if graph_genome.split('.')[-1]!='xg' and \
                graph_genome.split('.')[-1]!='vg':
                raise VGException("The genome graph must be in VG or XG format")
                die(1)
        
            elif graph_genome.split('.')[-1]=='vg': # we are given a vg genome
                vg=graph_genome
                xg=vgc.indexVG(vg)
            
            else: # we are given an xg genome
                xg=graph_genome
        
            without_vg_pipeline(cores, xg, bedfile, motif, bgfile, pseudocount, 
                                    pvalueT, no_reverse, dest, WITH_VG_CREATION) # add qvalue
            
        elif args.graph_genome_dir:
            
            gplus=True # defines if the input is a single xg or more than one
            
            without_vg_pipeline(cores, graph_genome_dir, bedfile, motif, bgfile, pseudocount,
                                    pvalueT, no_reverse, dest, WITH_VG_CREATION, gplus, chroms) # add qvalue
            
    else:
        # error in the pipeline flag
        msg="Wrong input. Unable to determine which pipeline to follow"
        raise PipelineException(msg)
        die(1)
        
    end=time.time()
    
    print('\nelapsed time', end-start)
    
    
### run main  point ###

if __name__=="__main__":
    main()
    
    
    