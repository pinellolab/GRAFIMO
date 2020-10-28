"""Functions to store GRAFIMO analysis results.

The results are stored in three different files:
* TSV file report, table of results, containing for each retrieved motif 
    occurrence candidates the analysis informations, like log-odds score,
    corresponding P-value and q-value, frequency in genome variation 
    graph haplotypes, etc.
* HTML file report, the same data stored in the TSV report are stored
    in an HTML file, which can be visualized with a commonly used
    web browser
* GFF file, can be loaded on the UCSC genome browser as custom track, to
    visualize on the gnome browser the motif candidates retrieved 
    through GRAFIMO analysis
"""

from grafimo.GRAFIMOException import VGException, SubprocessError, \
    NoDataFrameException, FileReadingException
from grafimo.utils import PHASE, TP, SOURCE, unique_lst, list_data, die, \
    DEFAULT_OUTDIR
from grafimo.workflow import Findmotif
from typing import List, Dict, Optional
import pandas as pd
import numpy as np
import subprocess
import time
import os


def write_results(results: pd.DataFrame,
                  motif_id: str,
                  motif_num: int,
                  args_obj: Findmotif
) -> None:
    """Write GRAFIMO analysis results in three files (TSV report, HTML 
    report, GFF3 file).

    To first two reports contain in a tabular format the motif occurrence
    candidates retrieved by GRAFIMO, with scores, P-values, q-values, etc.
    for each table entry.

    The third file contains a suitabl input for custom track on the UCSC
    genome browser.

    The user can also ask to display the results directly on the terminal
    enabling the correspondent flag, when calling GRAFIMO in command-line.

    Parameters
    ----------
    results : pandas.DataFrame
        results of GRAFIMO analysis
    motif_id : str
        motif ID
    motif_num : int
        number of searched motifs
    args_obj : Findmotif
        container for arguments needed to store the results
    """

    errmsg: str
    if not isinstance(args_obj, Findmotif):
        errmsg = "\n\nERROR: incorrect data-type. Exiting"
        raise ValueError(errmsg)

    if not isinstance(results, pd.DataFrame):
        errmsg = "\n\nERROR: results must be stored in a pandas DataFrame"
        raise NoDataFrameException(errmsg)

    # get resuls storing arguments
    outdir: str = args_obj.get_outdir()
    no_qvalue: bool = args_obj.get_no_qvalue()
    top_graphs: int = args_obj.get_top_graphs()
    verbose: bool = args_obj.get_verbose()

    vg: str
    if args_obj.has_graph_genome():
        vg = args_obj.get_graph_genome()
    elif args_obj.has_graph_genome_dir():
        vg = args_obj.get_graph_genome_dir()
    else:
        errmsg = "\n\nERROR: no VG given"
        raise VGException(errmsg)

    dirname_default: bool = False

    cwd: str = os.getcwd() 

    if outdir == DEFAULT_OUTDIR: 
        # to make unique the output directory we add the PID
        # to the name.
        #
        # This is useful when calling grafimo in different runs on the
        # same machine.

        # append the PID and the motif ID
        outdir = '_'.join(["grafimo_out", str(os.getpid()), motif_id])  
        dirname_default = True
    # end if

    cmd: str
    code: int
    if not os.path.isdir(outdir): 
        cmd = 'mkdir -p {0}'.format(outdir)
        code = subprocess.call(cmd, shell=True)

        if code != 0:
            errmsg = ''.join(["\n\nERROR: An error occurred while executing ", 
                              cmd, "\n"])
            raise SubprocessError(errmsg)

        os.chdir(outdir)

    else:
        os.chdir(outdir) 
        # NB the content of the directory will be automatically 
        # overwritten
    # end if

    print("\nWriting results in %s\n" % outdir)

    prefix: str
    # get the filename prefix
    if not dirname_default and motif_num > 1:
        # each file is labeled with the motif ID
        prefix = '_'.join(['grafimo_out', motif_id])
    else:
        prefix = 'grafimo_out'

    if verbose:
        start_tsv: float = time.time()

    # write the TSV
    results.to_csv(''.join([prefix, '.tsv']), sep='\t', encoding='utf-8')

    if verbose:
        end_tsv: float = time.time()
        print(
            "%s.tsv written in %.2fs" % (prefix, (end_tsv - start_tsv))
        )
        start_html: float = time.time()

    # write the HTML
    results.to_html(''.join([prefix, '.html']))

    if verbose:
        end_html: float = time.time()
        print(
            "%s.html written in %.2fs" % (prefix, (end_html - start_html))
        )
        start_gff: float = time.time()

    # write the GFF3
    writeGFF3(prefix, results, no_qvalue)

    if verbose:
        end_gff: float = time.time()
        print(
            "%s.gff written in %.2fs" % (prefix, (end_gff - start_gff))
        )

    # get the graphs of the top n regions
    regions: List[str] 
    if top_graphs > 0:
        regions = results['sequence_name'].to_list()
    
        # the results are empty
        if len(regions) == 0:
            print(
                "WARNING: no region was available. Are your results empty?\n"
            )
            os.chdir(cwd)
            return

        # get the n different top regions' graph
        regions = unique_lst(regions, size=top_graphs) 

        if verbose:
            print("Extracting %d region variation graphs" % top_graphs)

        if len(regions) < top_graphs:
            top_graphs = len(regions)
            print(
                "WARNING: possible to visualize only the top %d regions\n" % 
                (top_graphs)
                )

        # create the directory for the regions pictures
        image_dir: str
        if motif_num > 1:
            image_dir = '_'.join(["top_graphs", motif_id])
        else:
            image_dir = "top_graphs"

        if verbose:
            print("Graphs will be stored in %s" % image_dir)

        cmd = "mkdir -p {0}".format(image_dir)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = ' '.join(["\n\nERROR: an error occurred while executing",
                               cmd, "\n"])
            raise SubprocessError()

        os.chdir(image_dir)

        print("Writing the top %d graphs in %s\n" % (top_graphs, image_dir))

        for i in range(top_graphs):
            region: str = regions[i]
            # the VG accepts only queries like 1:100-200 and not like
            # chr1:100-200
            region = region.split("chr")[1]

            if verbose:
                print("Computing the PNG image of %s" % region)

            getRegion_graph(region, vg)

    os.chdir(cwd)

# end of write_results()


def writeGFF3(prefix: str, data: pd.DataFrame, no_qvalue: bool) -> None:
    """Write a GFF3 file 
    (https://www.ensembl.org/info/website/upload/gff3.html) containing 
    all the motif occurrence candidates retrieved by GRAFIMO.
    
    The resulting file can be loaded on the UCSC genome browser to view 
    the occurrence
        
    Parameters
    ----------
    prefix : str
        filename prefix
    data : pandas.DataFrame
        results of GRAFIMO analysis
    no_qvalue : bool
        if set to True the GFF3 entrues will contain also the 
        corresponding q-value
    """

    qvalue: bool
    errmsg: str
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

        data_list: List[List[str], List[str], List[str], List[int], List[int], 
                    List[str], List[float], List[float], List[str], List[int], 
                    List[str], Optional[List[float]]] 
        data_list = list_data(data, qvalue)

        if qvalue and len(data_list) < 12:
            errmsg = "\n\nERROR: wrong data size. Unable to write the GFF3 report\n"
            raise Exception(errmsg)

        data_list_size: int = len(data_list[0])
        for i in range(data_list_size):

            seqname: str = data_list[2][i]
            chrom: str = seqname.split(':')[0]  # takes only the chromosome name
            score: float = round(data_list[6][i], 1)
            strand: str = data_list[5][i]

            start: int
            end: int
            if strand == '-':
                # keep forward strand coordinates
                start = data_list[4][i]
                end = data_list[3][i]
            else:
                start = data_list[3][i]
                end = data_list[4][i]

            motifID: str = data_list[0][i]
            motifName: str = data_list[1][i]
            pvalue: float = np.format_float_scientific(data_list[7][i], 
                                                       exp_digits=2)
            sequence: str = data_list[8][i]
            reference: str = data_list[10][i]

            if qvalue:
                qvalue: float = np.format_float_scientific(data_list[11][i], 
                                                           exp_digits=2)

            att1: str = ''.join(['Name=', motifID, '_', seqname, strand, ':', 
                            reference])
            att2: str = ''.join(["Alias=", motifName])
            att3: str = ''.join(["ID=", motifID, '-', motifName, '-', seqname])
            att4: str = ''.join(['pvalue=', str(pvalue)])
            att5: str = ''.join(['sequence=', sequence, ';\n'])

            atts: str
            if qvalue:
                attqv: str = ''.join(['qvalue=', str(qvalue)])
                atts = ';'.join([att1, att2, att3, att4, attqv, att5])
            else:
                atts = ';'.join([att1, att2, att3, att4, att5])

            gffline: str = '\t'.join([chrom, SOURCE, TP, start, end, str(score),
                                 strand, PHASE, atts])

            f.write(gffline)

        # end for
    except:
        errmsg = ''.join(["\n\nERROR: unable to open or write data on ", prefix, 
                          ".gff"])
        raise FileReadingException(errmsg)
    finally:
        f.close()  
    # end try
    
# end of writeGFF3()
 

def getRegion_graph(region: str, genome_loc: str) -> None:
    """Get a PNG image representing the queried region of the genome
    variation graph.

    Parameters
    ----------
    region : str
        region for which the PNG image will be obtained
    genome_loc : str
        path to the genome variation graph(s)
    """

    errmsg: str
    cmd: str
    code: int
    vg_region: str

    if genome_loc.split('.')[-1] == 'xg':  
        # we have the whole genome graph

        if not os.path.isfile(genome_loc):
            errmsg = ' '.join(["\n\nERROR: Unable to locate", genome_loc, "\n"])
            raise Exception(errmsg)

        # extract the region PNG image 
        vg_region = ''.join([".", region, ".vg"])
        cmd = "vg find -x {0} -E -p {1} > {2}".format(genome_loc, region, 
                                                      vg_region)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessError("\n\nERROR: unable to execute %s" % cmd)

    else:  
        # we have a directory containing the genome graphs

        if not os.path.isdir(genome_loc):
            raise Exception("\n\nERROR: nable to locate %s" % genome_loc)

        # if given separate VGs (as built by GRAFIMO) they are called 
        # chr1.xg, chr20.xg, chrX.xg, etc.
        xg: str = ''.join(["chr", region.split(':')[0], '.xg'])

        xg = os.path.join(genome_loc, xg)
    
        # extract the PNG image of the region
        vg_region = ''.join([".", region, ".vg"])
        cmd = "vg find -x {0} -E -p {1} > {2}".format(xg, region, vg_region)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessError("\n\nERROR: unable to execute %s" % cmd)
    # end if
        
    dot_region: str = ''.join([".", region, ".dot"])
    cmd = "vg view {0} -dp > {1}".format(vg_region, dot_region)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise SubprocessError("\n\nERROR: unable to execute %s" % cmd)

    png_image: str = ''.join([region, '.png'])
    cmd = "dot -Tpng {0} -o {1}".format(dot_region, png_image)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise SubprocessError("\n\nERROR: unable to execute %s" % cmd)

    # clean the directory from unused files
    cmd = "rm .*.vg"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise SubprocessError("\n\nERROR: unable to execute %s" % cmd)

    cmd = "rm .*.dot"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise SubprocessError("\n\nERROR: unable to execute %s" % cmd)

# end of getRegion_graph()


def print_results(results):
    """Print GRAFIMO results on terminal without storing them on 
    the three files (TSV, HTML, GFF3)

    Parameters
    ----------
    results : pandas.DataFrame
        GRAFIMO results
    """

    if not isinstance(results, pd.DataFrame):
        errmsg: str = "\n\nERROR: the results must be stored in a pandas DataFrame"
        raise NoDataFrameException(errmsg)

    # little hack in pd df parameters to avoid the weird default
    # print of a DataFrame (cut the majority of lines)
    pd.set_option("display.max_rows", len(results))
    print()  # newline
    print(results)
    pd.reset_option("display.max_rows")

# end of print_results()

