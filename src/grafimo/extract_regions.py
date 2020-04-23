"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Extraction of the regions defined in the input BED file from the
queried genome graphs.

Are conceptually extracted subgraphs, from the genome graph,
corresponding to the regions defined in the input BED file.

"""


from grafimo.GRAFIMOException import SubprocessError, NotValidFFException, FileReadingException, VGException
from grafimo.utils import die, correct_path, CHROMS_LIST, printProgressBar, sigint_handler
from grafimo.score_sequences import ResultTmp
from grafimo.motif import Motif
import multiprocessing as mp
import subprocess
import signal
import warnings
import tempfile
import time
import sys
import os


verbose = False  # global variable


def get_regions(motif,
                args_obj):
    """
        Compute all sequences of length L (L is the
        motif width) from the VG(s).
        The sequences are extracted from the regions defined
        in the input BED file.
        ----
        Parameters:
            motif (Motif) : motif to search on the VG
            args_obj (Findmotif) : object storing the arguments
                                    required to extract the
                                    regions defined in the BED
                                    file, from the VG(s)
        ----
        Return:
            sequence_loc (str) : location of the tmp files,
                                    containing the extracted
                                    sequences
    """

    # check the input arguments
    if not isinstance(motif, Motif):
        errmsg = "\n\nERROR: unknown motif object type"
        raise ValueError(errmsg)

    if args_obj.has_graph_genome():
        vg = args_obj.get_graph_genome()

        if not isGraph_genome_xg(vg):
            errmsg = "\n\nERROR: the genome variation graph is not in XG format"
            raise VGException(errmsg)
        # end if

    elif args_obj.has_graph_genome_dir():
        vg = args_obj.get_graph_genome_dir()

    else:
        raise VGException("\n\nERROR: the genome variation graph is missing")
    # end if

    bedfile = args_obj.get_bedfile()
    motif_width = motif.getWidth()
    chroms = args_obj.get_chroms()
    cores = args_obj.get_cores()

    global verbose
    verbose = args_obj.get_verbose()

    print("\nExtracting regions defined in", bedfile, "\n")

    # read the regions where search the motif occurrences from the given BED file
    regions = getBEDregions(bedfile)

    if verbose:
        print("\nFound", len(regions), "regions in", bedfile)

    if chroms:
        # user defined subset of the chromosomes
        chr_list = [''.join(['chr', c]) for c in chroms]
    else:
        # all the chromosomes
        chr_list = [''.join(['chr', c]) for c in CHROMS_LIST]
    # end if

    # create a tmp working directory
    tmpwd = tempfile.mkdtemp(prefix='grafimo_')

    # if the tmp directory name already exists remove it
    # this shouldn't happen, but to be sure
    if os.path.isdir(tmpwd):
        cmd = 'rm -rf {0}'.format(tmpwd)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            raise SubprocessError(' '.join(["an error occurred executing", cmd, ". Exiting"]))
    # end if

    cmd = 'mkdir -p {0}'.format(tmpwd)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise SubprocessError(' '.join(["an error occurred executing", cmd, ". Exiting"]))

    # get the new location of graphs wrt the tmp dir
    cwd = os.getcwd()

    # enter the tmp dir where store the extracted sequences
    os.chdir(tmpwd)

    if verbose:
        start_re = time.time()

    # redefine default SIGINT handler
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool = mp.Pool(processes=cores)  # use #cores processes
    signal.signal(signal.SIGINT, original_sigint_handler)  # overwrite the default SIGINT handler to exit gracefully
    # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python

    if args_obj.has_graph_genome_dir():

        # vg -> directory containing a set of VGs
        if vg[-1] == "/":
            pass
        else:
            vg = ''.join([vg, "/"])
        # end if

        queries = []  # set of queries

        for region in regions:
            chrom = region['chr']
            start = region['start']
            stop = region['stop']

            if chrom in chr_list:

                # the chromosome is among the ones to query
                region_index = ''.join([chrom, ':', str(start), '-', str(stop)])
                region_name = ''.join([chrom, '_', str(start), '-', str(stop)])
                seqs = correct_path('./', region_name, '.tsv')

                xg = ''.join([vg, chrom, '.xg'])

                if not os.path.exists(xg):
                    errmsg = ''.join(["\n\nERROR: unable to use ", xg, ". Exiting"])
                    raise FileNotFoundError(errmsg)

                query = 'vg find -x {0} -E -p {1} -K {2} > {3}'.format(xg, region_index, motif_width, seqs)
                queries.append(query)

        # extract regions
        try:

            # query the VGs
            res = (pool.map_async(get_seqs, queries))

            if not verbose:
                it = 0
                while (True):
                    if res.ready():
                        # when finished call for the last time printProgressBar()
                        printProgressBar(tot, tot, prefix='Progress:',
                                         suffix='Complete', length=50)
                        break
                    # end if
                    if it == 0:
                        tot = res._number_left

                    remaining = res._number_left
                    printProgressBar((tot - remaining), tot, prefix='Progress:',
                                     suffix='Complete', length=50)
                    time.sleep(2)
                    it += 1
                # end while
            # end if

            ret = res.get(60 * 60 * 60)  # does not ignore signals

        except KeyboardInterrupt:
            pool.terminate()
            sigint_handler()

        else:
            pool.close()

            if verbose:
                end_re = time.time()
                msg = ''.join(["Extracted all regions from VGs stored in ", vg, ", in ",
                               str(end_re - start_re), "s"])
                print(msg)
            # end if
        # end try

    elif args_obj.has_graph_genome():

        queries = []  # set of queries

        for region in regions:
            chrom = region['chr']
            start = region['start']
            stop = region['stop']

            if chrom in chr_list:

                # the chromosome is among the ones to query
                region_index = ''.join([chrom, ':', str(start), '-', str(stop)])
                region_name = ''.join([chrom, '_', str(start), '-', str(stop)])
                seqs = correct_path('./', region_name, '.tsv')

                if not os.path.exists(vg):
                    errmsg = ''.join(["\n\nERROR: unable to use ", vg, ". Exiting"])
                    raise FileNotFoundError(errmsg)

                query = 'vg find -x {0} -E -p {1} -K {2} > {3}'.format(vg, region_index, motif_width, seqs)
                queries.append(query)

        # extract regions
        try:

            # query the VGs
            res = (pool.map_async(get_seqs, queries))

            if not verbose:
                it = 0
                while (True):
                    if res.ready():
                        # when finished call for the last time printProgressBar()
                        printProgressBar(tot, tot, prefix='Progress:',
                                         suffix='Complete', length=50)
                        break
                    # end if
                    if it == 0:
                        tot = res._number_left

                    remaining = res._number_left
                    printProgressBar((tot - remaining), tot, prefix='Progress:',
                                     suffix='Complete', length=50)
                    time.sleep(2)
                    it += 1
                # end while
            # end if

            ret = res.get(60 * 60 * 60)  # does not ignore signals

        except KeyboardInterrupt:
            pool.terminate()
            sigint_handler()

        else:
            pool.close()

            if verbose:
                end_re = time.time()
                msg = ''.join(["Extracted all regions from VGs stored in ", vg, ", in ",
                               str(end_re - start_re), "s"])
                print(msg)
            # end if
        # end try

    else:
        raise Exception("\n\nERROR: do not know how to proceed". Exiting)
    # end if

    sequence_loc = os.getcwd()  # the extracted sequences are store in the cwd
    os.chdir(cwd)  # get back to the origin

    return sequence_loc

# end of get_regions()


def get_xg_loc(toXGpath):
    """
        Get the path to the directory that contains the whole genome XG
        ----
        Parameters:
            toXGpath(str) : path to the xg
        ----
        Returns:
            toXGpath (str) : path to the directory containing the xg
    """

    for i in range(1, len(toXGpath)):
        if toXGpath[-i] == '/':
            bp = -i
            break
    # end for

    toXGpath = toXGpath[:bp]

    return toXGpath
# end of get_xg_loc()


def get_seqs(query):
    """
        Get the sequences of length L inside the queried
        region (L = motif width)
        ----
        Prameters:
            query (str) : call to VG to extract all sequences of length L
                            (L = motif width) from the given VG
        ----
        Returns:
            None
    """

    if verbose:
        region = query.split('-p')[1].split('-K')[0]
        print("Extracting sequences from region:", region)

    code = subprocess.call(query, shell=True)  # perform query

    if verbose:
        if code != 0:
            warnmsg = ''.join(["A problem occurred during sequences extraction in ", region])
            warnings.warn(warnmsg)
        # end if
    # end if

# end of get_seqs()


def isGraph_genome_xg(vg):
    """
        Check if the given genome variation graph is in xg format
        ----
        Parameters:
            vg (str) : genome variation graph
        ----
        Returns:
            (bool)
    """

    if not isinstance(vg, str):
        raise VGException("\n\nERROR: Invalid path to the genome graph. Cannot proceed")

    if vg.split('.')[-1] == 'xg':
        return True
    elif vg.split('.') == 'vg':
        return False
    else:
        errmsg = "\n\nERROR: do not know what to do with the given genome graph. Only XG or VG format allowed"
        raise VGException(errmsg)

# end of isGraph_genome_xg()


def getBEDregions(bedfile):
    """
        Read the BED file with regions to analyze
        ----
        Parameters:
            bedfile (str) : path to the BED file
        ----
        Returns:
            regions (list) : regions defined in the BED file
    """

    if bedfile.split('.')[-1] != 'bed':  # not a BED file
        raise NotValidFFException("The given BED file is not in BED format")

    regions = []

    # start reading
    try:
        with open(bedfile, mode='r') as inbed:  # open the BED file in read only mode
            for line in inbed:
                chrom, start, stop = line.split('\t')[0:3]
                region = {'chr': chrom, 'start': start, 'stop': stop}
                regions.append(region)
            # end for
        # end with

    except:  # not able to read the BED file
        msg = ' '.join(["\n\nError: unable to read", bedfile])
        raise FileReadingException(msg)

    else:
        return regions

    finally:
        inbed.close()  # close the file stream
    # end try
# end of getBEDregions()

