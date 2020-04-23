"""

workflows classes
"""

from argparse import Namespace


class Workflow(object):
    """
        Workflow base class
    """

    def __init__(self):
        pass


class BuildVG(Workflow):
    """
        Buildvg workflow class. This class stores all
        the arguments needed for 'grafimo buildvg'
        worflow execution.
    """

    reference_genome = None
    vcf = None
    chroms = None
    cores = None
    verbose = False  # default option
    test = False  # used to test the buildvg worflow (manually set)

    def __init__(self, args):
        if not isinstance(args, Namespace):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.linear_genome, str):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.vcf, str):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.chroms, list):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.cores, int):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.out, str):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.verbose, bool):
            raise ValueError("Incorrect command line arguments object type")

        self.reference_genome = args.linear_genome
        self.vcf = args.vcf
        self.chroms = args.chroms
        self.cores = args.cores
        self.outdir = args.out
        self.verbose = args.verbose
    # end of __init__()

    def get_reference_genome(self):
        if self.reference_genome:
            return self.reference_genome
        else:
            raise ValueError("No reference genome found")

    def get_vcf(self):
        if self.vcf:
            return self.vcf
        else:
            raise ValueError("No VCF file found")

    def get_chroms(self):
        if self.chroms:
            return self.chroms
        else:
            raise ValueError("No chromosomes list found")

    def get_cores(self):
        if self.cores:
            return self.cores
        else:
            raise ValueError("Unknown number of cores to use")

    def get_outdir(self):
        if self.outdir:
            return self.outdir
        else:
            raise ValueError("Unknown output directory")

    def get_verbose(self):
        return self.verbose

    def get_test(self):
        return self.test
# end of BuildVG


class Findmotif(Workflow):
    """
        Findmotif workflow class. This class stores all
        the arguments needed for 'grafimo findmotif'
        workflow execution
    """

    graph_genome = None  # if True graph_genome_dir must be None
    graph_genome_dir = None  # if True graph_genome must be None
    bedfile = None
    motif = None
    chroms = None
    bgfile = None
    pseudo = None
    thresh = None
    outdir = None
    cores = None
    top_graphs = 0  # default
    no_qval = False  # default
    no_rev = False  # default
    text_only = False  # default
    qvalT = False  # default
    verbose = False  # default
    test = True  # used to test the buildvg worflow (manually set)

    # the following variables cannot be both True at the same time
    _has_graph_genome = False
    _has_graph_genome_dir = False

    def __init__(self, args):
        if not isinstance(args, Namespace):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.graph_genome, str):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.graph_genome_dir, str):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.bedfile, str):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.chroms, list):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.bgfile, str):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.pseudo, float):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.threshold, float):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.out, str):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.cores, int):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.top_graphs, int):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.no_qvalue, bool):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.no_reverse, bool):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.text_only, bool):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.qval_t, bool):
            raise ValueError("Incorrect command line arguments object type")

        if not isinstance(args.verbose, bool):
            raise ValueError("Incorrect command line arguments object type")

        if args.graph_genome:
            self.graph_genome = args.graph_genome
            self._has_graph_genome = True
            self._has_graph_genome_dir = False
        elif args.graph_genome_dir:
            self.graph_genome_dir = args.graph_genome_dir
            self._has_graph_genome_dir = True
            self._has_graph_genome = False

        if self._has_graph_genome and self._has_graph_genome_dir:
            errmsg = "\n\nERROR: cannot be given both a path to a VG and a directory containing them"
            raise PermissionError(errmsg)

        self.bedfile = args.bedfile
        self.motif = args.motif
        self.chroms = args.chroms
        self.bgfile = args.bgfile
        self.pseudo = args.pseudo
        self.thresh = args.threshold
        self.cores = args.cores
        self.outdir = args.out
        self.top_graphs = args.top_graphs
        self.no_qval = args.no_qvalue
        self.no_rev = args.no_reverse
        self.text_only = args.text_only
        self.qvalT = args.qval_t
        self.verbose = args.verbose

    # end of __init__()

    def get_graph_genome(self):
        if not self._has_graph_genome:
            raise ValueError("No genome graph available")
        else:
            return self.graph_genome

    def get_graph_genome_dir(self):
        if not self._has_graph_genome_dir:
            raise ValueError("No genome graph directory given")
        else:
            return self.graph_genome_dir

    def get_bedfile(self):
        if not self.bedfile:
            raise ValueError("No BED file found")
        else:
            return self.bedfile

    def get_motif(self):
        if not self.motif:
            raise ValueError("NO motif file (MEME or JASPAR format) found")
        else:
            return self.motif

    def get_chroms(self):
        if not self.chroms:
            raise ValueError("No chromosomes list found")
        else:
            return self.chroms

    def get_bgfile(self):
        if not self.bgfile:
            raise ValueError("No background file found")
        else:
            return self.bgfile

    def get_pseudo(self):
        if not self.pseudo:
            raise ValueError("No pseudocount found")
        else:
            return self.pseudo

    def get_threshold(self):
        if not self.thresh:
            raise ValueError("No threshold found")
        else:
            return self.thresh

    def get_outdir(self):
        if not self.outdir:
            raise ValueError("No output directory found")
        else:
            return self.outdir

    def get_cores(self):
        if not self.cores:
            raise ValueError("Unknown number of cores to use")
        else:
            return self.cores

    def get_top_graphs(self):
        if self.top_graphs < 0:
            raise ValueError("The number of graph images to retrieve cannot be negative")
        else:
            return self.top_graphs

    def get_no_qvalue(self):
        return self.no_qval

    def get_no_reverse(self):
        return self.no_rev

    def get_text_only(self):
        return self.text_only

    def get_qvalueT(self):
        return self.qvalT

    def get_verbose(self):
        return self.verbose

    def get_test(self):
        return self.test

    def has_graph_genome(self):
        return self._has_graph_genome

    def has_graph_genome_dir(self):
        return self._has_graph_genome_dir

    def set_xg(self, vg):
        if vg.split('.')[-1] != 'xg':
            raise ValueError("Only XG accepted")
        self.graph_genome = vg

# end of Findmotif
