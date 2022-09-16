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

from grafimo.grafimo_errors import VGError, SubprocessError, FileWriteError
from grafimo.utils import (
    dftolist, 
    die,
    exception_handler,
    PHASE, 
    TP, 
    SOURCE, 
    DEFAULT_OUTDIR 
)
from grafimo.workflow import Findmotif
from grafimo.motif import Motif

from typing import List, Dict, Optional
from colorama import Fore, init

import pandas as pd
import numpy as np

import subprocess
import time
import sys
import os


def write_results(
    results: pd.DataFrame,
    motif: Motif,
    motif_num: int,
    args_obj: Findmotif,
    debug: bool
) -> None:
    """Write three reports storing GRAFIMo analysis results (TSV report, 
    HTML report, GFF3 report).

    The TSV and HTML reports stores the potential motif occurrence found 
    by GRAFIMO in tables.

    The GFF3 report stores annotations for the found motif occurrence 
    candidates. The GFF3 report can be loaded to the UCSC genome browser
    to visualize the locations of the found motif occurrences.

    ...

    Parameters
    ----------
    results : pandas.DataFrame
        Results
    motif_id : Motif
        Motif
    motif_num : int
        Number of searched motifs
    args_obj : Findmotif
        Commandline arguments container 
    debug : bool
        Trace the full error stack

    Returns
    -------
    None
    """

    if not isinstance(results, pd.DataFrame):
        errmsg = f"Expected {type(pd.DataFrame).__name__}, got {type(results).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if len(results) == 0:
        errmsg = "No potential motif occurrence retreived.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(motif, Motif):
        errmsg = f"Expected {type(Motif.__name__)}, got {type(motif).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(motif_num, int):
        errmsg = f"Expected {int.__name__}, got {type(motif_num).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if motif_num <= 0:
        errmsg = "No motif searched. Probably something went wrong.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(args_obj, Findmotif):
        errmsg = f"Expected {type(Findmotif).__name__}, got {type(args_obj).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    # get resuls storing arguments
    outdir = args_obj.outdir
    no_qvalue = args_obj.noqvalue
    top_graphs = args_obj.top_graphs
    verbose = args_obj.verbose
    if args_obj.has_graphgenome(): 
        vg = args_obj.graph_genome
    elif args_obj.has_graphgenome_dir: 
        vg = args_obj.graph_genome_dir
    else:
        errmsg = "No genome variation graph given.\n"
        exception_handler(VGError, errmsg, debug)
    dirname_default = False
    cwd = os.getcwd() 
    if outdir == DEFAULT_OUTDIR: 
        # to make unique the output directory we add the PID
        # to the name.
        #
        # This is useful when calling grafimo in different runs on the
        # same machine.
        outdir = "_".join(["grafimo_out", str(os.getpid()), motif.motif_id])  
        dirname_default = True
    if not os.path.isdir(outdir): 
        cmd = f"mkdir -p {outdir}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = f"An error occurred while executing \"{cmd}\".\n"
            exception_handler(SubprocessError, errmsg, debug)
        os.chdir(outdir)
    else:
        os.chdir(outdir)  # overwrite the content of the directory
    print(f"\nWriting results in {outdir}.\n")
    if not dirname_default and motif_num > 1:
        # each file is labeled with the motif ID
        prefix = "_".join(["grafimo_out", motif.motif_id])  
    else:
        prefix = "grafimo_out"
    if verbose:
        start_tsv = time.time()
    # write the TSV
    results.to_csv(".".join([prefix, "tsv"]), sep="\t", encoding="utf-8")
    if verbose:
        end_tsv = time.time()
        print("%s.tsv written in %.2fs" % (prefix, (end_tsv - start_tsv)))
        start_html = time.time()
    # write the HTML
    results.to_html(".".join([prefix, "html"]))
    if verbose:
        end_html = time.time()
        print("%s.html written in %.2fs" % (prefix, (end_html - start_html)))
        start_gff = time.time()
    # write the GFF3
    writeGFF3(prefix, results, no_qvalue, debug)
    if verbose:
        end_gff = time.time()
        print("%s.gff written in %.2fs" % (prefix, (end_gff - start_gff)))
    # get the graphs of the top n regions
    if top_graphs > 0:
        regions = set()
        for r in results["sequence_name"].tolist():
            if len(regions) >= top_graphs: break  # abort loop
            regions.add(r)  # avoid repeated regions
        # regions = set(results["sequence_name"].tolist()[:top_graphs])
        if len(regions) == 0:
            errmsg = "No region obtained, the results seems to be empty.\n"
            exception_handler(ValueError, errmsg, debug)
        if len(regions) < top_graphs:
            warnmsg = f"WARNING: requested {top_graphs} regions, obtained {len(regions)}.\n"
            sys.stderr.write(Fore.YELLOW + warnmsg + Fore.RESET)
        if verbose:
            print(f"Extracting {len(regions)} region variation graphs")
        # create the directory for the regions images
        if motif_num > 1:
            image_dir = "_".join(["top_graphs", motif.motifID])
        else:
            image_dir = "top_graphs"
        if verbose: 
            print(f"Graphs will be stored in {image_dir}.")
        cmd = f"mkdir -p {image_dir}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = f"An error ocurred while executing {cmd}."
            exception_handler(SubprocessError, errmsg, debug)
        assert os.path.isdir(image_dir)
        os.chdir(image_dir)
        print(f"Writing the top {len(regions)} graphs in {image_dir}\n")
        try:
            for r in regions:
                if verbose: 
                    print(f"Computing the PNG image of {r}")
                if args_obj.has_graphgenome():
                    get_region_graph(
                        r, 
                        args_obj.chroms_prefix, 
                        args_obj.namemap,
                        debug, 
                        graph_genome=args_obj.graph_genome
                    )
                elif args_obj.has_graphgenome_dir():
                    get_region_graph(
                        r, 
                        args_obj.chroms_prefix, 
                        args_obj.namemap,
                        debug, 
                        graph_genome_dir=args_obj.graph_genome_dir
                    )
                else:
                    errmsg = "Unknown VG type. Unable to print regions PNG images.\n"
                    exception_handler(ValueError, errmsg, debug)
        except:
            errmsg = f"An error occurred while computing PNG image of {r}.\n"
            exception_handler(VGError, errmsg, debug)
    os.chdir(cwd)

# end of write_results()


def writeGFF3(
    prefix: str, data: pd.DataFrame, no_qvalue: bool, debug: bool
) -> None:
    """Write GFF3 file (https://www.ensembl.org/info/website/upload/gff3.html). 
    
    The GFF3 file annotates the potential motif occurrences found by 
    GRAFIMO. The report can be loaded as custom track to the UCSC genome 
    browser for results visualization.

    ...
        
    Parameters
    ----------
    prefix : str
        Out filename prefix
    data : pandas.DataFrame
        Results
    no_qvalue : bool
        Ignore q-values
    debug : bool
        Trace the full error stack

    Returns
    -------
    None
    """

    if not isinstance(prefix, str):
        errmsg = f"Expected {str.__name__}, got {type(prefix).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(data, pd.DataFrame):
        errmsg = f"Expected {type(pd.DataFrame).__name__}, got {type(data).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(no_qvalue, bool):
        errmsg = f"Expected {bool.__name__}, got {type(no_qvalue).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    data_list = dftolist(data, no_qvalue, debug)
    try:
        gfffn = ".".join([prefix, "gff"])
        ofstream = open(gfffn, mode="w+")
        header = "##gff-version 3\n"
        ofstream.write(header)
        if not no_qvalue and len(data_list) != 12:
            errmsg = "Q-values columns seems to be missing.\n"
            exception_handler(ValueError, errmsg, debug)
        if no_qvalue and len(data_list) != 11:
            errmsg = "Too many or too few columns.\n"
            exception_handler(ValueError, errmsg, debug)
        data_list_size = len(data_list[0])
        for i in range(data_list_size):
            seqname = data_list[2][i]
            chrom = seqname.split(":")[0]  # take only chromosome name
            score = round(data_list[6][i], 1)
            strand = data_list[5][i]
            if strand == "-":  # keep forward strand coordinates
                start = str(data_list[4][i])
                stop = str(data_list[3][i])
            else:
                start = str(data_list[3][i])
                stop = str(data_list[4][i])
            motifID = data_list[0][i]
            motifName = data_list[1][i]
            pvalue = np.format_float_scientific(data_list[7][i], exp_digits=2)
            sequence = data_list[8][i]
            reference = data_list[10][i]
            if not no_qvalue: 
                qvalue = np.format_float_scientific(
                    data_list[11][i], exp_digits=2
                )
            # gff line attributes
            att1 = "".join(
                ["Name=", motifID, "_", seqname, strand, ":", reference]
            )
            att2 = "=".join(["Alias", motifName])
            att3 = "=".join(["ID", motifID, "-", motifName, "-", seqname])
            att4 = "=".join(["pvalue=", str(pvalue)])
            att5 = "=".join(["sequence=", sequence, ";\n"])  # end of gff line
            if not no_qvalue:
                attqv = "=".join(["qvalue", str(qvalue)])
                atts = ";".join([att1, att2, att3, att4, attqv, att5])
            else: atts = ";".join([att1, att2, att3, att4, att5])
            # full gff line
            gffline = "\t".join(
                [chrom, SOURCE, TP, start, stop, str(score), strand, PHASE, atts]
            )
            ofstream.write(gffline)
    except:
        errmsg = f"An error ocurred while writing {gfffn}.\n"
        exception_handler(FileWriteError, errmsg, debug)
    finally:
        ofstream.close()  
    
# end of writeGFF3()
 

def get_region_graph(
    region: str, 
    chroms_prefix: str,
    namemap: dict,
    debug: bool, 
    graph_genome: Optional[str] = None,
    graph_genome_dir: Optional[str] = None,
) -> None:
    """Compute the PNG image of genomic regions encoded in genome variation 
    graph(s).

    ...

    Parameters
    ----------
    region : str
        Genomic region
    chroms_prefix : str
        Chromosome prefix
    namemap : dict
        Chromosome names map
    debug : bool
        Trace the full error stack
    graph_genome : str
        Path to genome variation graph
    graph_genome_dir : str
        Path to directory of genome variation graphs
    """

    if not isinstance(region, str):
        errmsg = f"Expected {str.__name__}, got {type(region).__name__}.\n"
        exception_handler(TypeError, errmsg.format(type(region).__name__), debug)
    if not isinstance(chroms_prefix, str):
        errmsg = f"Expected {str.__name__}, got {type(chroms_prefix).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(namemap, dict):
        errmsg = f"Expected {dict.__name__}, got {type(namemap).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if graph_genome is None and graph_genome_dir is None:
        errmsg = "graph_genome and graph_genome_dir cannot be both None.\n"
        exception_handler(ValueError, errmsg, debug)
    if graph_genome is not None and graph_genome_dir is not None:
        errmsg = "graph_genome and graph_genome_dir cannot be both not None.\n"
        exception_handler(ValueError, errmsg, debug)
    if graph_genome is not None:
        if not isinstance(graph_genome, str):
            errmsg = f"Expected {str.__name__}, got {type(graph_genome).__name__}.\n"
            exception_handler(TypeError, errmsg, debug)
        if not os.path.isfile(graph_genome):
            errmsg = f"Unable to locate {graph_genome}.\n"
            exception_handler(FileNotFoundError, errmsg, debug)
    if graph_genome_dir is not None:
        if not isinstance(graph_genome_dir, str):
            errmsg = f"Expected {str.__name__}, got {type(graph_genome_dir).__name__}.\n"
            exception_handler(TypeError, errmsg, debug)
        if not os.path.isdir(graph_genome_dir):
            errmsg = f"Unable to locate {graph_genome_dir}."
            exception_handler(FileNotFoundError, errmsg, debug)
    if graph_genome and graph_genome_dir is None: 
        has_graphgenome = True
    else: 
        has_graphgenome = False  # graph_genome is None and graph_genome_dir == True
    if has_graphgenome:  # single genome variation graph
        vgregion = "".join([".", region, ".vg"])
        cmd = f"vg find -x {graph_genome} -E -p {region} > {vgregion}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = f"An error occurred while executing {cmd}.\n"
            exception_handler(SubprocessError, errmsg, debug)
    else:  # has_graphgenome == False
        chrom = region.split(":")[0]
        if bool(namemap): chrom = namemap[chrom]
        chrname = "".join([chroms_prefix, chrom])
        xg = os.path.join(graph_genome_dir, ".".join([chrname, "xg"]))
        vgregion = "".join([".", region, ".vg"])
        cmd = f"vg find -x {xg} -E -p {region} > {vgregion}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = f"An error occurred while executing {cmd}.\n"
            exception_handler(SubprocessError, errmsg, debug)
    dotregion = "".join([".", region, ".dot"])
    cmd = f"vg view {vgregion} -dp > {dotregion}"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        errmsg = f"An error occurred while executing {cmd}.\n"
        exception_handler(SubprocessError, errmsg, debug)
    pngimage = ".".join([region, "png"])
    cmd = f"dot -Tpng {dotregion} -o {pngimage}"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        errmsg = f"An error occurred while executing {cmd}.\n"
        exception_handler(SubprocessError, errmsg, debug)
    # remove unused files
    cmd = "rm -rf .*.vg"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        errmsg = f"An error occurred while executing {cmd}.\n"
        exception_handler(SubprocessError, errmsg, debug)
    cmd = "rm -rf .*.dot"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        errmsg = f"An error occurred while executing {cmd}.\n"
        exception_handler(SubprocessError, errmsg, debug)

# end of get_region_graph()


def print_results(results: pd.DataFrame, debug: bool) -> None:
    """Print GRAFIMO results to stdout. It is printed the tab-separated result
    summary.

    ...

    Parameters
    ----------
    results : pandas.DataFrame
        Results
    debug : bool
        Trace the full error stack
    """

    if not isinstance(results, pd.DataFrame):
        errmsg = f"Expected {pd.DataFrame}, got {type(results).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    # little hack in pd df parameters to avoid the weird default
    # print of a DataFrame (cut the majority of lines)
    pd.set_option("display.max_columns", None)  # ensure that all columns are displayed
    print()  # newline
    print(results)
    pd.reset_option("display.max_rows")

# end of print_results()

