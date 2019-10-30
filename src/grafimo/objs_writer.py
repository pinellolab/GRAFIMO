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
        sys.exit(1)
        
    try:
        f=open("grafimo_out.gff", 'w+')
    
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
    