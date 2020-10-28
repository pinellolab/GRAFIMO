"""Workflow classes definiton.

GRAFIMO provides the user two main workfloes to follow:
* genome variation graph construction
* genome variation graph scan for occurrences of a given DNA motif

For each one of the two workflows a class is created, which carries
all the arguments needed while executing the two operations.
"""

from typing import List
from argparse import Namespace


class Workflow(object):
    """
    This class is the Workflow base class which is the parent class of 
    BuildVG and Findmotif classes.

    ...

    Attributes
    ----------
    None

    Methods
    -------
    None
    """

    def __init__(self):
        pass

# end of Workflow


class BuildVG(Workflow):
    """
    This class represent the buildvg GRAFIMO workflow.

    The class carries all the arguments needed while building a genome 
    variation graph with user data, through GRAFIMO.

    ...

    Attributes
    ----------
    _reference_genome : str
        reference genometo enrich with the given genomic variants
    _vcf : str
        VCF file containing the genomic variants
    _reindex : bool
        if True the VCF file will be indexed with tabix, even if an 
        index already exists
    _chroms : list
        chromosomes for which the VG will be built
    _chroms_num : int
        number of chromosomes
    _cores : int
        number of cores to use during VG construction
    _outdir : str
        output directory
    _verbose : bool
        if True will be printed iformation about the execution
    _test : bool
        manually set test flag value

    Methods
    -------
    get_reference_genome()
        return the reference genome
    get_vcf()
        return the VCF file
    get_reindex()
        return if the VCF file has to be indexed 
    get_chroms()
        get the list of the chromosomes to consider
    get_chroms_num()
        get the number of chromosomes to consider
    get_cores()
        get the number of cores to use during genome variation graph 
        construction
    get_outdir()
        get the output directory
    get_verbose()
        get verbose flag
    get_test()
        get test flag
    """

    #-------------------------------------------------------------------
    # BuilVG attributes
    #-------------------------------------------------------------------
    _reference_genome: str
    _vcf: str
    _reindex: bool
    _chroms: List
    _chroms_num: int
    _cores: int
    _outdir: str
    _verbose: bool
    _test: bool = False


    #-------------------------------------------------------------------
    # BuildVG methods
    #-------------------------------------------------------------------
    def __init__(self, args: Namespace):
        errmsg: str = "\n\nERROR: incorrect command line arguments object type"
        if not isinstance(args, Namespace):
            raise ValueError(errmsg)

        if not isinstance(args.linear_genome, str):
            raise ValueError(errmsg)

        if not isinstance(args.vcf, str):
            raise ValueError(errmsg)

        if not isinstance(args.reindex, bool):
            raise ValueError(errmsg)

        if not isinstance(args.chroms, list):
            raise ValueError(errmsg)

        if not isinstance(args.cores, int):
            raise ValueError(errmsg)

        if not isinstance(args.out, str):
            raise ValueError(errmsg)

        if not isinstance(args.verbose, bool):
            raise ValueError(errmsg)

        self._reference_genome = args.linear_genome
        self._vcf = args.vcf
        self._reindex = args.reindex
        self._chroms = args.chroms
        self._chroms_num = len(args.chroms)
        self._cores = args.cores
        self._outdir = args.out
        self._verbose = args.verbose
    

    def get_reference_genome(self) -> str:
        if self._reference_genome:
            return self._reference_genome
        else:
            raise ValueError("No reference genome found")


    def get_vcf(self) -> str:
        if self._vcf:
            return self._vcf
        else:
            raise ValueError("No VCF file found")


    def get_reindex(self) -> bool:
        return self._reindex


    def get_chroms(self) -> List:
        if self._chroms is None:
            raise ValueError("No chromosomes list found")
        else:
            return self._chroms


    def get_chroms_num(self) -> int:
        if self._chroms_num == -1:
            raise ValueError("Forbidden number of chromosomes found")
        else:
            return self._chroms_num


    def get_cores(self) -> int:
        if self._cores:
            return self._cores
        else:
            raise ValueError("Unknown number of cores to use")

    def get_outdir(self) -> str:
        if self._outdir:
            return self._outdir
        else:
            raise ValueError("Unknown output directory")


    def get_verbose(self) -> bool:
        return self._verbose


    def get_test(self) -> bool:
        return self._test

# end of BuildVG


class Findmotif(Workflow):
    """
    This class represent the findmotif GRAFIMO workflow.

    The class carries all the arguments needed while scanning a genome 
    variation graph for the occurrences of aDNA motif, in the given
    set of genomic coordinates.

    ...

    Attributes
    ----------
    _graph_genome : str
        whole genome variation graph
    _graph_genome_dir : str
        directory containing a set of VGs 
    _bedfile : str
        BED file with the genomic coordinates
    _motif : str
        DNA motif in MEME or JASPAR format
    _chroms : list
        chromosomes to scan for motif occurrences
    _chroms_num : int
        number of chromosomes to scan for motif occurrences
    _bgfile : str
        background probabilities distribution file
    _pseudo : float
        pseudocount to add to motif matrix values
    _thresh : float
        threshold to filter reported motif occurrences
    _outdir : str
        output directory
    _cores : int
        cores to use during GRAFIMO scan
    _recomb : bool
        if True will be reported also those sequences obtained with the 
        given set of genomic variants, but never observed in the samples
        haplotypes
    _top_graphs : int
        number of regions to display as PNG images
    _no_qvalue : bool
        if True q-values will not be computed
    _no_rev : bool
        if True will be reported only sequences beloning to the forward
        strand
    _text_only : bool
        if True the results will be displayed on the screen and not 
        stored in files
    _qvalT : bool
        if True the threshold will be applied on q-values rather on 
        P-values
    _verbose : bool
        print additional information
    _test : bool
        test flag value
    _has_graph_genome : bool
        if True the Findmotif instance has a whole genome variation 
        graph and _has_graph_genome_dir must be False
    _has_graph_genome_dir : bool]
        if True the Findmotif instance has adirectory conatining genome 
        variation graphs and _has_graph_genome must be False

    Methods
    -------
    get_graph_genome():
        return the whole genome graph
    get_graph_genome_dir()
        return the directory containing genome variation graphs]
    get_befile()
        return the BED file
    get_motif()
        return the motif
    get_chroms()
        return the chromosomes to scan for motif occurrences
    get_chroms_num()
        return the number of chromosomes to scan for motif occurrences
    get_bgfile()
        return the background probability distribution file
    get_pseudo()
        return the pseudo count
    get_thresh()
        return the threshold to filter reported entries
    get_outdir()
        return the output directory
    get_cores()
        return the cores to use
    get_recomb()
        return recomb flag value
    get_top_graphs()
        return the number of regions to visualize as PNG images
    get_no_qvalue()
        return the no_qvalue flag value
    get_no_reverse()
        return the no_reverse flag value
    get_text_only()
        retur the text_only flag value
    get_qvalueT()
        return the qvalueT flag value
    get_verbose()
        return the verbose flag value
    get_test()
        return the test flag value
    has_graph_genome()
        return True if the current Findmotif instance has a whole genome 
        variation graph
    has_graph_genome_dir()
        return True if the current Findmotif instance has a directory
        conatining VGs
    set_xg(vg: str)
        set the whole genome variation graph 
    """

    #-------------------------------------------------------------------
    # Findmotif atrtributes
    #-------------------------------------------------------------------
    _graph_genome: str = None
    _graph_genome_dir: str = None
    _bedfile: str
    _motif: str
    _chroms: List
    _chroms_num: int = -1
    _bgfile: str
    _pseudo: float
    _thresh: float
    _outdir: str
    _cores: int
    _recomb = None
    _top_graphs: int
    _no_qvalue: bool
    _no_rev: bool
    _text_only: bool
    _qvalueT: bool
    _verbose: bool
    _test: bool = False  # used to test the findmotif worflow (manually set)

    # the following variables cannot be both True at the same time
    _has_graph_genome: bool = False
    _has_graph_genome_dir: bool = False


    #-------------------------------------------------------------------
    # Findmotif methods
    #-------------------------------------------------------------------
    def __init__(self, args):
        errmsg: str = "\n\nERROR: incorrect command line arguments object type"
        if not isinstance(args, Namespace):
            raise ValueError(errmsg)

        if not isinstance(args.graph_genome, str):
            raise ValueError(errmsg)

        if not isinstance(args.graph_genome_dir, str):
            raise ValueError(errmsg)

        if not isinstance(args.bedfile, str):
            raise ValueError(errmsg)

        if not isinstance(args.chroms, list):
            raise ValueError(errmsg)

        if not isinstance(args.bgfile, str):
            raise ValueError(errmsg)

        if not isinstance(args.pseudo, float):
            raise ValueError(errmsg)

        if not isinstance(args.threshold, float):
            raise ValueError(errmsg)

        if not isinstance(args.out, str):
            raise ValueError(errmsg)

        if not isinstance(args.cores, int):
            raise ValueError(errmsg)

        if not isinstance(args.recomb, bool):
            raise ValueError(errmsg)

        if not isinstance(args.top_graphs, int):
            raise ValueError(errmsg)

        if not isinstance(args.no_qvalue, bool):
            raise ValueError(errmsg)

        if not isinstance(args.no_reverse, bool):
            raise ValueError(errmsg)

        if not isinstance(args.text_only, bool):
            raise ValueError(errmsg)

        if not isinstance(args.qval_t, bool):
            raise ValueError(errmsg)

        if not isinstance(args.verbose, bool):
            raise ValueError(errmsg)

        if args.graph_genome:
            self._graph_genome = args.graph_genome
            self._has_graph_genome = True
            self._has_graph_genome_dir = False
        elif args.graph_genome_dir:
            self._graph_genome_dir = args.graph_genome_dir
            self._has_graph_genome_dir = True
            self._has_graph_genome = False

        if self._has_graph_genome and self._has_graph_genome_dir:
            errmsg = "\n\nERROR: cannot be given both a path to a VG and a "
            errmsg += "directory containing them"
            raise PermissionError(errmsg)

        self._bedfile = args.bedfile
        self._motif = args.motif
        self._chroms = args.chroms
        self._chroms_num = len(args.chroms)
        self._bgfile = args.bgfile
        self._pseudo = args.pseudo
        self._thresh = args.threshold
        self._cores = args.cores
        self._outdir = args.out
        self._recomb = args.recomb
        self._top_graphs = args.top_graphs
        self._no_qvalue = args.no_qvalue
        self._no_rev = args.no_reverse
        self._text_only = args.text_only
        self._qvalueT = args.qval_t
        self._verbose = args.verbose


    def get_graph_genome(self) -> str:
        if not self._has_graph_genome:
            raise ValueError("No genome graph available")
        else:
            return self._graph_genome


    def get_graph_genome_dir(self) -> str:
        if not self._has_graph_genome_dir:
            raise ValueError("No genome graph directory given")
        else:
            return self._graph_genome_dir


    def get_bedfile(self) -> str:
        if not self._bedfile:
            raise ValueError("No BED file found")
        else:
            return self._bedfile


    def get_motif(self) -> str:
        if not self._motif:
            raise ValueError("NO motif file (MEME or JASPAR format) found")
        else:
            return self._motif


    def get_chroms(self) -> List:
        if self._chroms is None:
            raise ValueError("No chromosomes list found")
        else:
            return self._chroms


    def get_chroms_num(self) -> int:
        if self._chroms_num == -1:
            raise ValueError("Forbidden number of chromosomes found")
        else:
            return self._chroms_num


    def get_bgfile(self) -> str:
        if not self._bgfile:
            raise ValueError("No background file found")
        else:
            return self._bgfile


    def get_pseudo(self) -> float:
        if not self._pseudo:
            raise ValueError("No pseudocount found")
        else:
            return self._pseudo


    def get_threshold(self) -> float:
        if not self._thresh:
            raise ValueError("No threshold found")
        else:
            return self._thresh


    def get_outdir(self) -> str:
        if not self._outdir:
            raise ValueError("No output directory found")
        else:
            return self._outdir


    def get_cores(self) -> int:
        if not self._cores:
            raise ValueError("Unknown number of cores to use")
        else:
            return self._cores


    def get_recomb(self) -> bool:
        return self._recomb


    def get_top_graphs(self) -> int:
        if self._top_graphs < 0:
            raise ValueError("The number of graph images to retrieve cannot be negative")
        else:
            return self._top_graphs


    def get_no_qvalue(self) -> bool:
        return self._no_qvalue


    def get_no_reverse(self) -> bool:
        return self._no_rev


    def get_text_only(self) -> bool:
        return self._text_only


    def get_qvalueT(self) -> bool:
        return self._qvalueT


    def get_verbose(self) -> bool:
        return self._verbose


    def get_test(self) -> bool:
        return self._test


    def has_graph_genome(self) -> bool:
        return self._has_graph_genome


    def has_graph_genome_dir(self) -> bool:
        return self._has_graph_genome_dir


    def set_xg(self, vg: str) -> None:
        if vg.split('.')[-1] != 'xg':
            raise ValueError("Only XG accepted")
        self._graph_genome = vg

# end of Findmotif

