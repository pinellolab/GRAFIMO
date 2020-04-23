"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Construction and indexing of the genome-graphs.
The reference used is given by the user as like as the VCF
used to infer the variants and the alternative paths on
the graphs.

We choose to split the genome graph for the genome in a graph for
each chromosome, in order to reduce the amount of memory and
computational resources needed.

NB There are NO drawbacks using this approach instead of creating the whole genome graph

"""


from grafimo.utils import CHROMS_LIST, die
from grafimo.workflow import BuildVG
from grafimo.GRAFIMOException import VGException, ValueException, SubprocessError
import subprocess
import time
import os


def get_reference_genome_from_ucsc():
    """
        Download the reference genome (hg38 assembly), from the UCSC
        database, in the current directory and returns the path to it
        ----
        Parameters:
            None
        ----
        Returns:
            genome (str) : path to the genome downloaded (in .fa format)
    """

    # download genome
    address = 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
    cmd = 'wget -c {0}'.format(address)
    code = subprocess.call(cmd, shell=True)  # downloaded in the current directory

    if code != 0:
        errmsg = ''.join(["\n\nERROR: an error occurred while executing ", cmd, ". Exiting"])
        raise SubprocessError(errmsg)
        die(1)

    # decompress genome
    print("Uncompressing the genome...")

    genome_comp = './hg38.fa.gz'

    cmd = 'gunzip {0}'.format(genome_comp)
    code = subprocess.call(cmd, shell=True)

    if code != 0:
        errmsg = ''.join(["\n\nERROR: an error occurred while executing ", cmd, ". Exiting"])
        raise SubprocessError(errmsg)
        die(1)

    # remove FASTA.GZ file if still present
    if os.path.exists(genome_comp):
        cmd = 'rm {0}'.format(genome_comp)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            errmsg = ''.join(["\n\nERROR: an error occurred while executing ", cmd, ". Exiting"])
            raise SubprocessError(errmsg)
            die(1)

    # get the path to the genome file
    genome_uncomp = "./hg38.fa"
    genome = os.path.abspath(genome_uncomp)

    return genome

# end of get_reference_genome_from_ucsc()


def get_1000GProject_vcf():
    """
        Downloads a VCF file from the 1000 Genome Project database,
        containing SNPs and indels for each subject involved in the
        study and returns the path to it
        ----
        Parameters:
            None
        ----
        Returns:
            vcf (str) : path to the vcf downloaded vcf file (in .vcf.gz)
    """

    # download the VCF
    address = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/'
    address += '1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/'
    address += 'ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz'

    cmd = 'wget -c {0}'.format(address)

    code = subprocess.call(cmd, shell=True)

    if code != 0:
        errmsg = ''.join(["\n\nERROR: An error occurred while executing ", cmd, ". Exiting"])
        raise SubprocessError(errmsg)
        die(1)

    vcf_file = './ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz'
    vcf = os.path.abspath(vcf_file)

    return vcf

# end of get1000GProject_vcf()


def construct_vg(buildvg_args):
    """
        Create the genome graph, for the reference and VCF file given
        in input by the user.

        The genome is not built as a single whole genome graph but a
        single graph is constructed for each chromosome.
        This choice was made to avoid memory issues and make able
        also the less powerful machines to run GRAFIMO.

        There is NO drawback using this approach wrt
        construct the whole genome graph and query it.
        ----
        Parameters:
            chroms (list) : list of chromosomes for whicgh the genome
                            graph will be constructed
            linear_genome (str) : path to the linear genome used as
                                    reference to build the genome
                                    graphs
            vcf (str) : path to the VCF file used to build the genome
                        graphs
        ----
        Return:
            None
    """

    if not isinstance(buildvg_args, BuildVG):
        raise ValueError("Unknown arguments object type. Cannot Build the genome variation graph. Exiting")
        die(1)

    # read the arguments to build the VGs
    chroms = buildvg_args.get_chroms()
    threads = buildvg_args.get_cores()
    outdir = buildvg_args.get_outdir()
    verbose = buildvg_args.get_verbose()
    test = buildvg_args.get_test()

    if test:
        reference = get_reference_genome_from_ucsc()
        vcf = get_1000GProject_vcf()

    else:
        reference = buildvg_args.get_reference_genome()
        vcf = buildvg_args.get_vcf()
    # end if

    if verbose:
        print("using reference genome: ", reference)
        print("Using VCF file: ", vcf, "\n\n")
    # end if

    cwd = os.getcwd()

    # check if the VCF file has already been indexed with tabix
    if not tbiexist(vcf):
        msg = ''.join(["TBI file not found for ", vcf.split('/')[-1], ". Indexing the VCF file with tabix..."])
        print(msg)
        cmd = 'tabix -p vcf {0}'.format(vcf)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            # tabix didn't work
            errmsg = ''.join(["\n\nERROR: an error occurred while executing ", cmd, ". Exiting"])
            raise SubprocessError(errmsg)
            die(1)

    else:  # update the indexed VCF
        msg = ''.join(["Reindexing ", vcf.split('/')[-1], "..."])
        print(msg)

        # remove the existing TBI file
        cmd = "rm {0}".format(''.join([vcf, ".tbi"]))
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            errmsg = ''.join(["\n\nERROR: an error occurred while executing ", cmd, ". Exiting"])
            raise SubprocessError(errmsg)
            die(1)

        # reindex the VCF
        cmd = "tabix -p vcf {0}".format(vcf)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            # tabix didn't work
            errmsg = ''.join(["\n\nERROR: an error occurred while executing ", cmd, ". Exiting"])
            raise SubprocessError(errmsg)
            die(1)
        # end if
    # end if

    # enter the output directory
    os.chdir(outdir)

    # build the VG for each chromosome or a user defined
    # subset of them
    for chrom_n in chroms:

        chrom = ''.join(['chr', chrom_n])  # to call vg construct we need both the
                                           # chromosome number and it preceded by 'chr

        vg = chrom + '.vg'

        # build the VG for the current chromosome
        if verbose:
            start_build = time.time()

        code = build_vg(vg, reference, vcf, chrom, chrom_n, threads)
        if code != 0:
            msg = '\n\nERROR: an error occurred during {0} construction. '.format(vg)
            msg += 'Unable to build the VG of the genome using {0} and {1}'.format(reference,
                                                                                   vcf)
            raise VGException(msg)
            die(1)
        # end if

        if verbose:
            end_build = time.time()
            msg = "Elapsed time to build {0} ".format(vg)
            msg = ''.join([msg, str(end_build - start_build), "s"])
            print(msg)
        # end if

        # to query efficiently the VGs we index them (VG -> XG)
        if verbose:
            start_index = time.time()

        msg = ''.join(["Indexing ", vg, '...'])
        print(msg)

        code = indexVG(vg, threads)

        if code != 0:
            errmsg = "\n\nERROR: an error occurred during indexing {0}.\nUnable to index {0}. Exiting".format(vg)
            raise VGException(errmsg)
            die(1)
        # end if

        if verbose:
            end_index = time.time()
            msg = "Elapsed time to index {0} ".format(vg)
            msg = ''.join([msg, str(end_index - start_index), "s"])
            print(msg)
        # end if

        # The majority of applications work only with indexed graph,
        # so to save disk space is worth to delete the VGs and keep
        # only the XGs (is simple to get back using VG built-in functions)
        if verbose:
            print("Deleting {0}".format(vg))

        cmd = 'rm {0}'.format(vg)
        subprocess.call(cmd, shell=True)

        if code != 0:  # we have errors in the vg indexing
            errmsg = ''.join(["\n\nERROR: an error occurred while executing ", cmd, ". Exiting"])
            raise SubprocessError()
            die(1)
        # end if
    # end for

    # get the VGs location
    graphs_loc = os.getcwd()

    # return to the original working directory
    os.chdir(cwd)

# end of construct_vg()


def build_vg(vg, linear_genome, vcf, chrom_name, chrom_num, threads):
    """
        Build the genome graph of the given chromosome
        ----
        Parameters:
            vg (str) : output name of the genome graph
            linear_genome (str) : path to the linear genome location
            vcf (str) : path to the VCF file
            chrom_name (str) : chromosome name as used by UCSC (e.g. chr1)
            chrom_num (int) : chromosome number
            threads (int) : number of threads to use during VG construction
        ----
        Returns:
            done (int) : if an error occurred during the building of
                            vg is returned 1, 0 otherwise
    """

    done = 0

    build = 'vg construct -C -R {0} -p -t {6} -n {1}={2} -r {3} -v {4} > {5}'.format(chrom_num, chrom_num,
                                                                                     chrom_name, linear_genome, vcf,
                                                                                     vg, threads)

    code = subprocess.call(build, shell=True)

    if code != 0:
        # an error occurred during the construction of the VG for
        # the current chromosome
        done = 1

    return done

# end of build_vg()


def indexVG(vg, threads):
    """
        Index the vg genome graph [vg => xg]. This step is required to
        efficiently extract sequences from VG, in defined regions.
        ----
        Parameters:
            vg (str) : path to the VG genome graph
            threads (int) : number of threads to use during VG indexing
        ----
        Returns:
            done (int) : if an error occurred during the indexing of the
                            given VG is returned 1, 0 otherwise
    """

    if not isinstance(vg, str):
        raise ValueError("Invalid path to the genome variation graph. Exiting")

    if not os.path.exists(vg):
        errmsg = "unable to find {0}".format(vg)
        raise FileExistsError(errmsg)

    done = 0

    # take chromosome name and add it the XG extension
    xg = vg.split('.')[-2]
    xg += '.xg'

    # perform indexing of the current genome variation graph
    vg_index = 'vg index -p -t {2} --progress -x {0} {1}'.format(xg, vg, threads)
    code = subprocess.call(vg_index, shell=True)

    if code != 0:
        done = 1

    return done

# end of indexVG()


def tbiexist(vcf):
    """
        Check if the given VCF file has already been indexed
        using tabix
        ----
        Parameters:
            vcf (str) : path to the VCF file, for which it will
                            be searched a TBI file in the
                            current directory
        ----
        Returns:
            (bool)
    """

    vcf_tbi = ''.join([vcf, '.tbi'])

    if os.path.isfile(vcf_tbi):
        return True

    return False

# end of tbiexist()

