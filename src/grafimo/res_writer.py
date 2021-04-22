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

from grafimo.GRAFIMOException import VGError, SubprocessError, FileWriteError
from grafimo.utils import PHASE, TP, SOURCE, dftolist, die, \
    DEFAULT_OUTDIR, exception_handler
from grafimo.workflow import Findmotif
from grafimo.motif import Motif

from typing import List, Dict, Optional

import pandas as pd
import numpy as np

import subprocess
import time
import os


def write_results(
    results: pd.DataFrame,
    motif: Motif,
    motif_num: int,
    args_obj: Findmotif,
    debug: bool
) -> None:
    """Write GRAFIMO results in three files (TSV report, HTML report, GFF3 file).

    The TSV and HTML reports stores the found potential motif occurrence in
    tabular format

    The GFF3 report stores annotations for the found motif occurrence candidates.

    ...

    Parameters
    ----------
    results : pandas.DataFrame
        analysis results
    motif_id : Motif
        motif
    motif_num : int
        number of searched motifs
    args_obj : Findmotif
        commandline arguments container 
    debug : bool
        trace the full error stack
    """

    if not isinstance(results, pd.DataFrame):
        errmsg = "Expected pandas.DataFrame, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(results).__name__), debug)
    if len(results) == 0:
        errmsg = "No potential motif occurrence retreived.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(motif, Motif):
        errmsg = "Expected Motif, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif).__name__), debug)
    if not isinstance(motif_num, int):
        errmsg = "Expected int, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_num).__name__), debug)
    if motif_num <= 0:
        errmsg = "No motif searched. Probably something went wrong.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(args_obj, Findmotif):
        errmsg = "Expected Findmotif, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(args_obj).__name__), debug)

    # get resuls storing arguments
    outdir: str = args_obj.outdir
    no_qvalue: bool = args_obj.noqvalue
    top_graphs: int = args_obj.top_graphs
    verbose: bool = args_obj.verbose
    if args_obj.has_graphgenome(): vg = args_obj.graph_genome
    elif args_obj.has_graphgenome_dir: vg = args_obj.graph_genome_dir
    else:
        errmsg = "No genome variation graph given.\n"
        exception_handler(VGError, errmsg, debug)
    dirname_default: bool = False
    cwd: str = os.getcwd() 
    if outdir == DEFAULT_OUTDIR: 
        # to make unique the output directory we add the PID
        # to the name.
        #
        # This is useful when calling grafimo in different runs on the
        # same machine.
        outdir = "_".join(["grafimo_out", str(os.getpid()), motif.motifID])  
        dirname_default = True
    if not os.path.isdir(outdir): 
        cmd = "mkdir -p {}".format(outdir)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = "An error occurred while executing {}.\n"
            exception_handler(SubprocessError, errmsg.format(cmd), debug)
        os.chdir(outdir)
    else:
        os.chdir(outdir)  # overwrite the content of the directory
    print("\nWriting results in %s.\n" % outdir)
    if not dirname_default and motif_num > 1:
        prefix = "_".join(["grafimo_out", motif.motifID])  # each file is labeled with the motif ID
    else:
        prefix = "grafimo_out"
    if verbose:
        start_tsv: float = time.time()
    # write the TSV
    results.to_csv(".".join([prefix, "tsv"]), sep='\t', encoding='utf-8')
    if verbose:
        end_tsv: float = time.time()
        print("%s.tsv written in %.2fs" % (prefix, (end_tsv - start_tsv)))
        start_html: float = time.time()
    # write the HTML
    results.to_html(".".join([prefix, "html"]))
    if verbose:
        end_html: float = time.time()
        print("%s.html written in %.2fs" % (prefix, (end_html - start_html)))
        start_gff: float = time.time()
    # write the GFF3
    writeGFF3(prefix, results, no_qvalue, debug)
    if verbose:
        end_gff: float = time.time()
        print("%s.gff written in %.2fs" % (prefix, (end_gff - start_gff)))
    # get the graphs of the top n regions
    if top_graphs > 0:
        regions = set(results["sequence_name"].tolist()[:top_graphs])
        if len(regions) == 0:
            errmsg = "No region obtained, the results seems to be empty.\n"
            exception_handler(ValueError, errmsg, debug)
        if len(regions) < top_graphs:
            warnmsg = "WARNING: requested %d regions, obtaned %d.\n"
            print(warnmsg % (top_graphs, len(regions)))
        if verbose:
            print("Extracting %d region variation graphs" % len(regions))
        # create the directory for the regions images
        if motif_num > 1:
            image_dir = "_".join(["top_graphs", motif.motifID])
        else:
            image_dir = "top_graphs"
        if verbose: print("Graphs will be stored in %s." % image_dir)
        cmd = "mkdir -p {0}".format(image_dir)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = "An error ocurred while executing {}."
            exception_handler(SubprocessError, errmsg.format(cmd), debug)
        assert os.path.isdir(image_dir)
        os.chdir(image_dir)
        print("Writing the top %d graphs in %s\n" % (len(regions), image_dir))
        try:
            for r in regions:
                if verbose: print("Computing the PNG image of {}".format(region))
                get_region_graph(region, vg)
        except:
            errmsg = "An error occurred while computing PNG image of {}.\n"
            exception_handler(VGError, errmsg.format(r), debug)
    os.chdir(cwd)

# end of write_results()


def writeGFF3(prefix: str, data: pd.DataFrame, no_qvalue: bool, debug: bool) -> None:
    """Write GFF3 file (https://www.ensembl.org/info/website/upload/gff3.html). 
    
    The GFF3 file annotates the potential motf occurrences found by GRAFIMO. The
    report can be loaded as custom track to the UCSC genome browser for results
    visualization.

    ...
        
    Parameters
    ----------
    prefix : str
        out filename prefix
    data : pandas.DataFrame
        analysis results
    no_qvalue : bool
        ignore q-values
    debug : bool
        trace the full error stack
    """

    if not isinstance(prefix, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format_map(type(prefix).__name__), debug)
    if not isinstance(data, pd.DataFrame):
        errmsg = "Expected pandas.DataFrame, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(data).__name__), debug)
    if not isinstance(no_qvalue, bool):
        errmsg = "Expected bool, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(no_qvalue).__name__), debug)

    data_list = dftolist(data, no_qvalue, debug)
    try:
        gfffn = ".".join([prefix, "gff"])
        ofstream = open(gfffn, mode='w+')
        header = "##gff-version 3\n"
        ofstream.write(header)
        if not no_qvalue and len(data_list) != 12:
            errmsg = "Q-values columns seems to be missing.\n"
            exception_handler(ValueError, errmsg, debug)
        if no_qvalue and len(data_list) != 11:
            errmsg = "Too many or too few columns.\n"
            exception_handler(ValueError, errmsg, debug)
        data_list_size: int = len(data_list[0])
        for i in range(data_list_size):
            seqname: str = data_list[2][i]
            chrom: str = seqname.split(":")[0]  # take only chromosome name
            score: float = round(data_list[6][i], 1)
            strand: str = data_list[5][i]
            if strand == "-":  # keep forward strand coordinates
                start = str(data_list[4][i])
                stop = str(data_list[3][i])
            else:
                start = str(data_list[3][i])
                stop = str(data_list[4][i])
            motifID: str = data_list[0][i]
            motifName: str = data_list[1][i]
            pvalue: float = np.format_float_scientific(data_list[7][i], exp_digits=2)
            sequence: str = data_list[8][i]
            reference: str = data_list[10][i]
            if not no_qvalue: 
                qvalue: float = np.format_float_scientific(data_list[11][i], exp_digits=2)
            # gff line attributes
            att1: str = "".join(
                ["Name=", motifID, "_", seqname, strand, ":", reference]
            )
            att2: str = "".join(["Alias=", motifName])
            att3: str = "".join(["ID=", motifID, "-", motifName, "-", seqname])
            att4: str = "".join(["pvalue=", str(pvalue)])
            att5: str = "".join(["sequence=", sequence, ";\n"])  # end of gff line
            if not no_qvalue:
                attqv: str = "".join(["qvalue=", str(qvalue)])
                atts = ";".join([att1, att2, att3, att4, attqv, att5])
            else: atts = ";".join([att1, att2, att3, att4, att5])
            # full gff line
            gffline: str = "\t".join(
                [chrom, SOURCE, TP, start, stop, str(score), strand, PHASE, atts]
            )
            ofstream.write(gffline)
    except:
        errmsg = "An error ocurred while writing {}.\n"
        exception_handler(FileWriteError, errmsg.format(gfffn), debug)
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
        genomic region
    chroms_prefix : str
        chromosome prefix
    namemap : dict
        chromosome names map
    debug : bool
        trace the full error stack
    graph_genome : str
        path to genome variation graph
    graph_genome_dir : str
        path to directory of genome variation graphs
    """

    if not isinstance(region, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(region).__name__), debug)
    if not isinstance(chroms_prefix, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(chroms_prefix).__name__), debug)
    if not isinstance(namemap, dict):
        errmsg = "Expected dict, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(namemap).__name__), debug)
    if graph_genome is None and graph_genome_dir is None:
        errmsg = "graph_genome and graph_genome_dir cannot be both None.\n"
        exception_handler(ValueError, errmsg, debug)
    if graph_genome is not None and graph_genome_dir is not None:
        errmsg = "graph_genome and graph_genome_dir cannot be both not None.\n"
        exception_handler(ValueError, errmsg, debug)
    if graph_genome is not None:
        if not isinstance(graph_genome, str):
            errmsg = "Expected str, got {}.\n"
            exception_handler(TypeError, errmsg.format(type(graph_genome).__name__), debug)
        if not os.path.isfile(graph_genome):
            errmsg = "Unable to locate {}.\n"
            exception_handler(FileNotFoundError, errmsg.format(graph_genome), debug)
    if graph_genome_dir is not None:
        if not isinstance(graph_genome_dir, str):
            errmsg = "Expected str, got {}.\n"
            exception_handler(TypeError, errmsg.format(type(graph_genome_dir).__name__), debug)
        if not os.path.isdir(graph_genome_dir):
            errmsg = "Unable to locate {}."
            exception_handler(FileNotFoundError, errmsg.format(graph_genome_dir), debug)

    if graph_genome and graph_genome_dir is None: has_graphgenome = True
    else: has_graphgenome = False  # graph_genome is None and graph_genome_dir == True

    if has_graphgenome:  # single genome variation graph
        vgregion = "".join([".", region, ".vg"])
        cmd = "vg find -x {} -E -p {} > {}".format(graph_genome, region, vgregion)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = "An error occurred while executing {}.\n"
            exception_handler(SubprocessError, errmsg.format(cmd), debug)
    else:  # has_graphgenome == False
        chrom = region.split(":")[0]
        if bool(namemap): chrom = namemap[chrom]
        chrname = "".join([chroms_prefix, chrom])
        xg = os.path.join(graph_genome_dir, ".".join([chrname, "xg"]))
        vgregion = "".join([".", region, ".vg"])
        cmd = "vg find -x {} -E -p {} > {}".format(xg, region, vgregion)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            errmsg = "An error occurred while executing {}.\n"
            exception_handler(SubprocessError, errmsg.format(cmd), debug)
    dotregion = "".join([".", region, ".dot"])
    cmd = "vg view {} -dp > {}".format(vgregion, dotregion)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        errmsg = "An error occurred while executing {}.\n"
        exception_handler(SubprocessError, errmsg.format(cmd), debug)
    pngimage = ".".join([region, "png"])
    cmd = "dor -Tpng {} -o {}".format(dotregion, pngimage)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        errmsg = "An error occurred while executing {}.\n"
        exception_handler(SubprocessError, errmsg.format(cmd), debug)
    # remove unused files
    cmd = "rm -rf .*.vg"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        errmsg = "An error occurred while executing {}.\n"
        exception_handler(SubprocessError, errmsg.format(cmd), debug)
    cmd = "rm -rf .*.dot"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        errmsg = "An error occurred while executing {}.\n"
        exception_handler(SubprocessError, errmsg.format(), debug)

# end of get_region_graph()


def print_results(results: pd.DataFrame, debug: bool):
    """Print GRAFIMO results to stdout. It is printed the tab-separated result
    summary.

    ...

    Parameters
    ----------
    results : pandas.DataFrame
        analysis results
    debug : bool
        trace the full error stack
    """

    if not isinstance(results, pd.DataFrame):
        errmsg = "Expected pandas.DataFrame, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(results).__name__), debug)

    # little hack in pd df parameters to avoid the weird default
    # print of a DataFrame (cut the majority of lines)
    pd.set_option("display.max_rows", len(results))
    print()  # newline
    print(results)
    pd.reset_option("display.max_rows")

# end of print_results()

