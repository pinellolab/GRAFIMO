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
from . import handle_exception as he

def writeresults(objs, dest):
    
    if not isinstance(objs, list):
        code=he.throw_not_list_err()
        sys.exit(code)
        
    if len(objs) <= 0:
        code=he.throw_objs_list_length()
        sys.exit(code)
        
    cwd=os.getcwd()
    
    if not os.path.isdir(dest): # tthe directory not exist
    
        cmd='mkdir {0}'.format(dest) # create the out directory
        code=subprocess.call(cmd, shell=True)
    
        if code!=0:
            raise Exception('error while executing a shell command (mkdir)')
            sys.exit(1)
            
        os.chdir(dest)
            
    else:
        os.chdir(dest) # for some reason it exists
        
        # the content will be automatically rewritten
    
    # write objects in dest
    
    for obj in objs:
        
        if isinstance(obj, pd.DataFrame):
            df=obj
            
            # write the tsv
            df.to_csv('grafimo_out.tsv', sep='\t', encoding='utf-8')
            
            # write the html
            df.to_html('grafimo_out.html')
            
            # write the gff
            writeGFF3(df)
            
            
        
        ### if matplotlib istance ###
        ## write plots in dest ##
        
    
    os.chdir(cwd)
    
def writeGFF3(data):
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
        raise Exception("DataFrame given is not an instance of pandas.DataFrame")
        sys.exit(1)
        
    try:
        f=open("grafimo_out.gff", mode='w+')
    
        header="##gff-version 3\n"
        f.write(header)
    
        data_idxs=list(data.index)
        
        for idx in data_idxs:
            
            seqid=data.loc[idx, 'sequence_name']
            source='grafimo'
            tp='nucleotide_motif'
            start=data.loc[idx, 'start']
            end=data.loc[idx, 'stop']
            score=round(data.loc[idx, 'score'], 1)
            strand=data.loc[idx, 'strand']
            phase='.'
        
            motifID=data.loc[idx, 'motif_id']
            motifName=data.loc[idx, 'motif_alt_id']
            pvalue=data.loc[idx, 'p-value']
            sequence=data.loc[idx, 'matched_sequence']
            att1=''.join(['Name=', motifID, '_', seqid, strand])
            att2=''.join(["Alias=", motifName])
            att3=''.join(["ID=", motifID, '-', motifName, '-', seqid])
            att4=''.join(['pvalue=', str(pvalue)])
            att5=''.join(['sequence=', sequence])
            atts=':'.join([att1, att2, att3, att4, att5])
        
            gffline='\t'.join([seqid, source, tp, start, end, str(score), 
                                   strand, phase, atts, ';\n'])
    
            f.write(gffline)
            
    except IOError:
        print('Error: unable to open or write data on grafimo_out.gff')
        sys.exit(1)
        
#def writeHTML(data):
#    
#    if not isinstance(data, pd.DataFrame):
#        raise Exception("DataFrame given is not an instance of pandas.DataFrame")
#        sys.exit(1)
#        
#    try:
#        f=open("grafimo_out.html", mode='w+')
#        
#        header=getHTMLHeader()
#        f.write(header)
#        f.write(grafimoHTMLBody())
#        
#    except IOError:
#        print('Error: unable to open or write data on grafimo_out.html')
#        sys.exit(1)
#        
#        
#def grafimo_version():
#
#    v='<p>\n'
#    v+=' '.join(["GRAFIMO version", __version__, "\n"])
#    v+='</p>\n'
#    
#    return v
#
#def getHTMLHeader():
#    
#    head='<head>\n'
#    head+='<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">\n'
#    head+='<meta charset="UTF-8">\n'
#    head+='<title>GRAFIMO Results Summary</title>\n'
#    head+='<style type="text/css">\n'
#    head+='td.left {text-align: left;}\n'
#    head+='td.right {text-align: right; padding-right: 1cm;}\n'
#    head+='</style>\n'
#    head+='</head>\n'
#
#    return head
#
#def grafimo_citation():
#    
#    citation='<p>\n'
#    citation+='For some informations on how interpret the results you can look at <a href="https://github.com/InfOmics/GRAFIMO/README.md">https://github.com/InfOmics/GRAFIMO/README.md</a>.<br>\n'
#    citation+='To get a copy of GRAFIMO you can download and build it from <a href="https://github.com/InfOmics/GRAFIMO">https://github.com/InfOmics/GRAFIMO</a>\n'
#    citation+='</p>\n'
#    citation+='<p>\n'
#    citation+='If you use GRAFIMO in your research, please cite:<br>\n'
#    citation+='XXX, "GRAFIMO", <i>xxx</i>, <b>x</b>(x):xxxx-xxxx, xxxx.\n'
#    citation+='<a href="http://xxxx">[full text]</a>\n'
#    citation+='</p>\n'
#
#    return citation
#
#def grafimoHTMLMotif(motif, bedfile):
#    
#    m='<hr>\n'
#    m+='<div style="padding-left: 0.75in; line-height: 1em; font-family: monospace;">\n'
#    m+=getHTMLBED(bedfile)
#    m+=getHTMLMotif
#<p>
#  MOTIFS /Users/manuel/Desktop/LM-MedicalBioinformatics/Thesis/2019-08-16/PWMs/MA0139.1.meme (DNA)
#  <table>
#    <thead>
#      <tr>
#        <th style="border-bottom: 1px dashed;">MOTIF</th>
#        <th style="border-bottom: 1px dashed; padding-left: 1em;">WIDTH</th>
#        <th style="border-bottom: 1px dashed; padding-left: 1em;text-align:left;" >
#         BEST POSSIBLE MATCH
#        </th>
#      </tr>
#    </thead>
#    <tbody>
#      <tr>
#        <td style="text-align:right;">MA0139.1</td>
#        <td style="text-align:right;padding-left: 1em;">19</td>
#        <td style="text-align:left;padding-left: 1em;">TGGCCACCAGGGGGCGCTA</td>
#       </tr>
#    </tbody>
#  </table>
#</p>
#<p>
#Random model letter frequencies (/Users/manuel/Desktop/bg_nt):
#<br/>
#
#A 0.295 C 0.205 G 0.205 T 0.295 </p>
#</div>
#
#def getHTMLBED(bed):
#    
#    b='<p>\n'
#    b+=' '.join(["BED FILE", bed, "\n"])
#    b+='</p>\n'
#    b+='<br />\n'
#    
#    BEDstats=computeBEDstats(bed)
#    
#    b+=' '.join(['Regions extracted:', BEDstats)
#    b+='</p>\n'
#    
#    return b
#
#def getHTMLMotif(motif):
#    
#    
#
#def grafimoHTMLBody(options):
#    
#    motif=options.getMotif()
#    bedfile=options.getBED()
#    
#    body='<body color="#000000">\n'
#    body+='<a name="top buttons"></a>\n'
#    body+='<hr>\n'
#    body+='<table summary="buttons" align="left" cellspacing="0">\n'
#    body+='<tr>\n'
#    body+='<td><a href="#motif_and_regions"><b>Motif and Regions</b></a></td>\n'
#    body+='<td><a href="#sec_1"><b>High-scoring Motif Occurences</b></a></td>\n'
#    body+='<td><a href="#submission_infos"><b>Submission Informations</b></a></td>\n'
#    body+='<td><a href="grafimo_out.tsv"><b>Results in TSV Format</b></a></td>\n'
#    body+='<td><a href="grafimo_out.gff"><b>Results in GFF3 Format</b></a></td>\n'
#    body+='</tr>\n'
#    body+='</table>\n'
#    body+='<br/>\n'
#    body+='<br/>\n'
#    body+='<hr/>\n'
#    body+='<center><big><b>GRAFIMO - GRAph-based Find Individual Occurrences</b></big></center>\n'
#    body+='<hr>\n'
#    body+=grafimo_version()
#    body+=grafimo_citation()
#    body+=grafimoHTMLMotif(motif, bedfile)
# 
#        
        