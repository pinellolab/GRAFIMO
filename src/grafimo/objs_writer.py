"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Scripts that writes the results to a user defined (or the default defined) 
directory.

"""

import sys
import os
import subprocess
import pandas as pd
import numpy as np
from grafimo.GRAFIMOException import ValueException, SubprocessException, NoDataFrameException, FileReadingException
from grafimo.utils import PHASE, TP, SOURCE, die, unique_lst, list_data

def writeresults(objs, dest, motifID, top_graphs, genome_loc, qvalue):
    
    if not isinstance(objs, list):
        raise ValueException("The results to write must be in a list")
        die(1)

    if len(objs) <= 0:
        raise Exception("Empty results list. Cannot proceed")
        die(1)
        
    cwd = os.getcwd()

    dest = '_'.join([dest, motifID])
    
    if not os.path.isdir(dest): # tthe directory not exist
    
        cmd = 'mkdir -p {0}'.format(dest) # create the out directory
        code = subprocess.call(cmd, shell=True)
    
        if code != 0:
            raise SubprocessException(' '.join(["An error occurred while executing", cmd,
                                                    ". Exiting"]))
            die(1)
            
        os.chdir(dest)
            
    else:
        os.chdir(dest) # for some unknown reason it exists
        
        # the content will be automatically rewritten
    
    # write objects in dest
    printWriteResultsMsg(dest)
    
    for obj in objs:
        
        if isinstance(obj, pd.DataFrame):
            df = obj
            
            # write the tsv
            df.to_csv('grafimo_out.tsv', sep='\t', encoding='utf-8')
            
            # write the html
            df.to_html('grafimo_out.html')
            
            # write the gff
            writeGFF3(df, qvalue)

    # get the graphs of the top n regions 
    if top_graphs > 0:
        # get the region list (it's already sorted)
        regions = df['sequence_name'].to_list()
        regions = unique_lst(regions, size=top_graphs)

        if len(regions) < top_graphs:
            top_graphs = len(regions)
            print("Warning: possible to visualize only the top " + str(top_graphs) + " regions")

        # create the directory for the images of the regions
        image_dir = "top_graphs"
        cmd = "mkdir -p {0}".format(image_dir)

        code = subprocess.call(cmd, shell = True)

        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

        # enter the newly created directory
        os.chdir(image_dir)

        print("Writing the top " + str(top_graphs) + " graphs in " + image_dir)
        print()

        for i in range(top_graphs):
            region = regions[i]

            getRegion_graph(region, genome_loc)

    os.chdir(cwd)

    
def writeGFF3(data, qvalue):
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
        raise NoDataFrameException("DataFrame given is not an instance of pandas.DataFrame")
        die(1)
        
    try:
        f = open("grafimo_out.gff", mode='w+')
    
        header = "##gff-version 3\n"
        f.write(header)

        data_list = list_data(data, qvalue)

        if qvalue and len(data_list) < 11:
            raise Exception("Some data were lost while transforming the result dataframe in a matrix")

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
        raise FileReadingException("Unable to open or write data on grafimo_out.gff")
        die(1)

    finally:
        f.close() # close the file stream


def getRegion_graph(region, genome_loc):
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

    if genome_loc.split('.')[-1] == 'xg': # we have the whole genome graph
        
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
        code = subprocess.call(cmd, shell = True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

        cmd = "rm *.dot"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise SubprocessException("Error while executing " + cmd)
            die(1)

    else: # we have a directory containing the genome graphs

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
        code = subprocess.call(cmd, shell = True)
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
        directly to the terminal
        ----
        Parameters:
            results (pd.DataFrame) : dataframe containing the analysis
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
    print() # newline
    print(results)
    pd.reset_option("display.max_rows")


def printWriteResultsMsg(dest):

    print()
    msg = ' '.join(["Writing results in", dest])
    print(msg)
    print()

