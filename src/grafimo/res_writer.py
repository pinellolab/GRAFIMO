"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Scripts that writes the results to a user defined (or the default defined)
directory.

"""

from grafimo.GRAFIMOException import VGException, SubprocessError, NoDataFrameException, FileReadingException
from grafimo.utils import PHASE, TP, SOURCE, unique_lst, list_data
from grafimo.workflow import Findmotif
import pandas as pd
import numpy as np
import subprocess
import time
import os


def write_results(results,
                  motif_id,
                  motif_num,
                  args_obj):

    if not isinstance(args_obj, Findmotif):
        errmsg = "\n\nERROR: incorrect data-type. Exiting"
        raise ValueError(errmsg)

    if not isinstance(results, pd.DataFrame):
        errmsg = "\n\nERROR: results must be stored in a pandas DataFrame"
        raise NoDataFrameException(errmsg)

    # read arguments
    outdir = args_obj.get_outdir()
    no_qvalue = args_obj.get_no_qvalue()
    top_graphs = args_obj.get_top_graphs()
    verbose = args_obj.get_verbose()

    if args_obj.has_graph_genome():
        vg = args_obj.get_graph_genome()
    elif args_obj.has_graph_genome_dir():
        vg = args_obj.get_graph_genome_dir()
    else:
        raise VGException("\n\nERROR: no VG given")

    dirname_default = False

    cwd = os.getcwd()  # get current directory

    if outdir == "grafimo_out":  # default option
        outdir = '_'.join([outdir, str(os.getpid()), motif_id])  # append the PID and the motif ID
        dirname_default = True
    # end if

    if not os.path.isdir(outdir):  # the directory does not exist

        cmd = 'mkdir -p {0}'.format(outdir)  # create the out directory
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            errmsg = ''.join(["An error occurred while executing ", cmd, ". Exiting"])
            raise SubprocessError(errmsg)
        # end if

        os.chdir(outdir)

    else:
        os.chdir(outdir)  # it already exists
        # NB the content will be automatically overwritten
    # end if

    # write results files in outdir
    print("\nWriting results in", outdir, "\n")

    # get the filename prefix
    if not dirname_default and motif_num > 1:
        # each file is labeled with the motif ID
        prefix = '_'.join(['grafimo_out', motif_id])
    else:
        prefix = 'grafimo_out'

    if verbose:
        start_tsv = time.time()

    # write the TSV
    results.to_csv(''.join([prefix, '.tsv']), sep='\t', encoding='utf-8')

    if verbose:
        end_tsv = time.time()
        msg = ''.join([prefix, ".tsv written in ", str(end_tsv - start_tsv), "s"])
        print(msg)
        start_html = time.time()

    # write the HTML
    results.to_html(''.join([prefix, '.html']))

    if verbose:
        end_html = time.time()
        msg = ''.join([prefix, ".html written in ", str(end_html - start_html), "s"])
        print(msg)
        start_gff = time.time()

    # write the GFF3
    writeGFF3(prefix, results, no_qvalue)

    if verbose:
        end_gff = time.time()
        msg = ''.join([prefix, ".gff written in ", str(end_gff - start_gff), "s"])
        print(msg)

    # get the graphs of the top n regions
    if top_graphs > 0:
        # get the region list (it's already sorted)
        regions = results['sequence_name'].to_list()
        regions = unique_lst(regions, size=top_graphs)

        if verbose:
            print("Extracting", top_graphs, "region variation graphs")

        if len(regions) < top_graphs:
            top_graphs = len(regions)
            print("Warning: possible to visualize only the top " + str(top_graphs) + " regions")

        # create the directory for the images of the regions
        if motif_num > 1:
            image_dir = ''.join(["top_graphs", motif_id])
        else:
            image_dir = "top_graphs"

        if verbose:
            print("Graphs will be stored in", image_dir)

        cmd = "mkdir -p {0}".format(image_dir)

        code = subprocess.call(cmd, shell=True)

        if code != 0:
            raise SubprocessError("\n\nERROR: an error occurred while executing " + cmd)

        # enter the newly created directory
        os.chdir(image_dir)

        print("Writing the top " + str(top_graphs) + " graphs in " + image_dir)
        print()

        for i in range(top_graphs):
            region = regions[i]

            if verbose:
                print("Computing the PNG image of", region)

            getRegion_graph(region, vg)

    os.chdir(cwd)


def writeGFF3(prefix,
              data,
              no_qvalue):
    """
        Write a GFF file for the hits found by grafimo.
        Are followed the conventions given at
        https://www.ensembl.org/info/website/upload/gff3.html
        ----
        Parameters:
            data (pandas DataFrame) : DataFrame containing the hits from
                                      GRAFIMO pipline
        ----
        Returns:
            None
    """

    if not isinstance(data, pd.DataFrame):
        errmsg = "\n\nERROR: the object is not an instance of pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if no_qvalue:
        qvalue = False
    else:
        qvalue = True

    try:
        f = open(''.join([prefix, '.gff']), mode='w+')

        header = "##gff-version 3\n"
        f.write(header)

        data_list = list_data(data, qvalue)

        if qvalue and len(data_list) < 11:
            raise Exception("Some data were lost while transforming the result data frame in a matrix")

        data_list_size = len(data_list[0])

        for i in range(data_list_size):

            seqname = data_list[2][i]
            chrom = seqname.split(':')[0]  # takes only the chromosome name
            score = round(data_list[6][i], 1)
            strand = data_list[5][i]

            if strand == '-':
                # keep forward strand coordinates
                start = data_list[4][i]
                end = data_list[3][i]
            else:
                start = data_list[3][i]
                end = data_list[4][i]

            motifID = data_list[0][i]
            motifName = data_list[1][i]
            pvalue = np.format_float_scientific(data_list[7][i], exp_digits=2)
            sequence = data_list[8][i]
            reference = data_list[9][i]

            if qvalue:
                qvalue = np.format_float_scientific(data_list[10][i], exp_digits=2)

            att1 = ''.join(['Name=', motifID, '_', seqname, strand, ':', reference])
            att2 = ''.join(["Alias=", motifName])
            att3 = ''.join(["ID=", motifID, '-', motifName, '-', seqname])
            att4 = ''.join(['pvalue=', str(pvalue)])
            att5 = ''.join(['sequence=', sequence, ';\n'])

            if qvalue:
                attqv = ''.join(['qvalue=', str(qvalue)])
                atts = ';'.join([att1, att2, att3, att4, attqv, att5])
            else:
                atts = ';'.join([att1, att2, att3, att4, att5])

            gffline = '\t'.join([chrom, SOURCE, TP, start, end, str(score),
                                 strand, PHASE, atts])

            f.write(gffline)

        # end for

    except:
        errmsg = ''.join(["\n\nERROR: unable to open or write data on ", prefix, ".gff"])
        raise FileReadingException(errmsg)

    finally:
        f.close()  # close the file stream


def getRegion_graph(region,
                    genome_loc):
    """
        Extract the queried region from the graph genome
        ----
        Parameters:
            region (str) : region to extract
            genome_loc (str) : path to the genome
        ----
        Returns:
            None
    """

    if genome_loc.split('.')[-1] == 'xg':  # we have the whole genome graph

        if not os.path.isfile(genome_loc):
            raise Exception("Unable to locate genome " + genome_loc)
            die(1)

        png_file = ''.join([region, '.png'])

        # extract the PNG image of the region
        vg_region = ''.join([region, ".vg"])
        cmd = "vg find -x {0} -E -p {1} > {2}".format(genome_loc, region, vg_region)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

        dot_region = ''.join([region, ".dot"])
        cmd = "vg view {0} -dp > {1}".format(vg_region, dot_region)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

        png_file = ''.join([region, '.png'])
        cmd = "dot -Tpng {0} -o {1}".format(dot_region, png_file)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

        # clean the directory from unused files
        cmd = "rm *.vg"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

        cmd = "rm *.dot"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

    else:  # we have a directory containing the genome graphs

        if not os.path.isdir(genome_loc):
            raise Exception("Unable to locate directory " + genome_loc)
            die(1)

        xg = ''.join([region.split(':')[0], '.xg'])
        if genome_loc[-1] == '/':
            xg = ''.join([genome_loc, xg])
        else:
            xg = '/'.join([genome_loc, xg])

        # extract the PNG image of the region
        vg_region = ''.join([region, ".vg"])
        cmd = "vg find -x {0} -E -p {1} > {2}".format(xg, region, vg_region)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

        dot_region = ''.join([region, ".dot"])
        cmd = "vg view {0} -dp > {1}".format(vg_region, dot_region)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

        png_file = ''.join([region, '.png'])
        cmd = "dot -Tpng {0} -o {1}".format(dot_region, png_file)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

        # clean the directory from unused files
        cmd = "rm *.vg"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

        cmd = "rm *.dot"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

    # end if else


def print_results(results):
    """
        Print the DataFrame containing the results
        directly to the screen (no file is created)
        ----
        Parameters:
            results (pd.DataFrame) : data frame containing the analysis
                                        results
        ----
        Returns:
            None
    """

    if not isinstance(results, pd.DataFrame):
        raise NoDataFrameException("The results must be stored in a pandas DataFrame")

    # little hack in pd df parameters to avoid the weird default
    # print of a DataFrame (cut the majority of lines)
    pd.set_option("display.max_rows", len(results))
    print()  # newline
    print(results)
    pd.reset_option("display.max_rows")

