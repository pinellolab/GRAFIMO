"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

The script defines the Motif class, that represents a motif and other functions
to manage its construction 

"""

import pandas as pd
import numpy as np
import sys
from . import handle_exception as he

class Motif(object):
    """
        Class to represent a motif.
        Are stored:
            - the motif matrix
            - the motif matrix scaled
            - the reverse matrix
            - the reverse matrix scaled
            - the pvalue matrix
            - the width of the motif
            - the minimum value in the motif matrix (both forward and reverse)
            - the transcription factor's name
            - the alphabet on which the matrix is built
    """
    
    motif_matrix=pd.DataFrame()
    motif_matrix_scaled=pd.DataFrame()
    motif_pval_matrix=pd.DataFrame()
    width=-1
    minValue=-1
    motifID='' # jaspar ID
    motifName='' # usual name
    alphabet=''
    
    def __init__(self, motif_matrix, width, motifID='', motifName=''):
        
        if motif_matrix.empty:
            code=he.throw_empty_motif_error()
            sys.exit(code)
            
        if not isinstance(motif_matrix, pd.DataFrame):
            code=he.throw_motif_matrix_not_pd_dataframe()
            sys.exit(code)
            
        if not isinstance(width, int) or width < 0:
            code=he.throw_incorrect_width_error()
            sys.exit(code)
            
        if not isinstance(motifID, str):
            code=he.throw_not_str_error()
            sys.exit(code)
        
        self.motif_matrix=motif_matrix
        self.width=width
        self.motifID=motifID
        self.motifName=motifName
    
    def setMotif_matrix(self, motif_matrix):
        
        if motif_matrix.empty:
            code=he.throw_empty_motif_error()
            sys.exit(code)
            
        elif not isinstance(motif_matrix, pd.DataFrame):
            code=he.throw_motif_matrix_not_pd_dataframe()
            sys.exit(code)
            
        else:
            self.motif_matrix=motif_matrix
            
    def setMotif_matrix_scaled(self, scaled_matrix):
        
        if scaled_matrix.empty:
            code=he.throw_empty_motif_error()
            sys.exit(code)
            
        elif not isinstance(scaled_matrix, pd.DataFrame):
            code=he.throw_motif_matrix_not_pd_dataframe()
            sys.exit(code)
            
        else:
            self.motif_matrix_scaled=scaled_matrix
            
    def setMotif_pval_matrix(self, pval_mat):
        
        self.motif_pval_matrix=pval_mat
              
    def setWidth(self, width):
        
        if not isinstance(width, int):
            code=he.throw_incorrect_width_error()
            sys.exit(code)
            
        else:
            self.width=width
        
    def setMotifID(self, motifID):
        
        if not isinstance(motifID, str):
            code=he.throw_not_str_error()
            sys.exit(code)
            
        elif not motifID:
            code=he.throw_empty_TF_name_error()
            sys.exit(code)
            
        else:
            self.motifID=motifID
            
    def setMotifName(self, motifName):
        
        if not isinstance(motifName, str):
            code=he.throw_not_str_error()
            sys.exit(code)
            
        elif not motifName:
            code=he.throw_empty_TF_name_error()
            sys.exit(code)
            
        else:
            self.motifName=motifName
            
    def setAlphabet(self, alphabet):
        
        if not isinstance(alphabet, str):
            code=he.throw_not_dna_alphabet_error()
            sys.exit(code)
            
        else:
            self.alphabet=alphabet
        
    def getMotif_matrix(self):
        
        return self.motif_matrix
    
    def getMotif_matrix_scaled(self):
        
        return self.motif_matrix_scaled
    
    def getMotif_pval_mat(self):
        
        return self.motif_pval_matrix
    
    def getWidth(self):
        
        return self.width
    
    def getMotifID(self):
        
        return self.motifID
    
    def getMotifName(self):
        
        return self.motifName
    
    def getMinValue(self):
        
        return self.minValue
    
    def getAlphabet(self):
        
        return self.alphabet
    
    def has_scaled_matrix(self):
        
        if not self.motif_matrix_scaled.empty:
            return True
        
        else:
            return False
        
    def compute_minValue(self):
        
        motif_matrix=self.getMotif_matrix()
        minValue=motif_matrix.min().sum()
        self.minValue=minValue
        
    

def isJaspar_ff(motif_file):
    """
        Check if a file is in .jaspar format
        ----
        Parameters:
            motif_file (str) : path to the motif file
        ----
        Returns:
            (bool)
    """
    
    if motif_file and isinstance(motif_file, str):
        ff=motif_file.split('.')[-1]
        
        if ff=='jaspar':
            return True
        else:
            return False
        
    else:
        return False # the motif file was not given as a path or the path is of length 0
            
def build_motif_pwm_jaspar(motif_file, bg_file, pseudocount):
    """
        Build the motif object, given a .jaspar motif file
        ----
        Parameters:
            motif_file (str) : path to the motif file
            bg_file (str) : path to the background file
            pseudocount (float) : value to add to the motif counts
        ----
        Returns:
            motif (Motif) : returns a Motif object
    """
    
    # check if the input file is in jaspar format
    if not isJaspar_ff(motif_file):
        code=he.throw_not_jaspar_error()
        sys.exit(code)
    
    # build the motif pwm
    nucs = []
    fqs = []
    bgs=readBGfile(bg_file)
    
    try:
        mf=open(motif_file, mode='r') # open the file in reading only mode
        
        header=str(mf.readline()[1:]) # read the header
        motifID, motifName=header.split('\t')[0:2] # get the jaspar ID and the usual
                                             # name for the motif
        motifName=motifName[:-1]
                                             
        for line in mf:
            line = line.strip()
            nuc = line.strip()[:1] # get the nucleotide
            fq = list(map(float, line.strip()[1:].split()[1:][:-1])) # get the frequency
            
            nucs.append(nuc)
            fqs.append(fq)
        
    except: # something went wrong
        code=he.throw_motif_file_not_read_error()
        sys.exit(code)
            
    finally:
        mf.close() # close the motif file
        
    alphabet=''.join(nucs)
            
    motif_count = pd.DataFrame(fqs, index=nucs) # raw counts
    motif_width=len(motif_count.columns)
    
    motif_prob = (motif_count/motif_count.sum(0)) # probabilities
    motif_log_odds=pd.DataFrame(index=nucs, columns=range(motif_width), data=0)
   
    for j in range(motif_width): 
        
        site_counts=sum(motif_count.loc[:,j])
        total_counts=site_counts+pseudocount
        for nuc in alphabet:
            bg=bgs[nuc]
            freq = ((pseudocount*bg)+(motif_prob.loc[nuc,j]*site_counts))/total_counts
            
            if freq <= 0:
                freq = 0.0000005
            if bg <= 0:
                bg = 0.0000005
            
            motif_log_odds.loc[nuc,j]=np.log2(freq/bg)
            
    motif=Motif(motif_log_odds, motif_width, motifID=motifID, 
                    motifName=motifName) # create the motif
    motif.setAlphabet(alphabet)
            
    return motif

def scale_pwm(motif_matrix):
    """
        Computes the scaled matrix, given a motif matrix
        ----
        Parameters:
            motif_matrix (str) : pat to the motif file
        ----
        Returns:
            motif_matrix (pd.DataFrame) : scaled motif matrix
    """
    
    if not isinstance(motif_matrix, pd.DataFrame):
        code=he.throw_motif_matrix_not_pd_dataframe()
        sys.exit(code)
        
    else:
    
        min_val=min(motif_matrix.min())
        motif_matrix=motif_matrix - min_val
    
        max_val=max(motif_matrix.max())
        scale_factor=100.0/ float(max_val)
        motif_matrix=round(motif_matrix * scale_factor)
    
        motif_matrix[:]=motif_matrix[:].astype(int)
     
        return motif_matrix

def comp_pval_mat(motif):
    """
        Computes the p-value matrix needed for the comoputation of the p-value,
        using the Dynamic Programming algorithm
        ----
        Parameters:
            motif (Motif) : motif object
        ----
        Returns:
            motif (Motif) : motif object updated with the p-value matrix
    """    
    if not isinstance(motif, Motif):
        code=he.throw_not_Motif_instance()
        sys.exit(code)
        
    elif not motif.has_scaled_matrix():
        code=he.throw_scaled_matrix_not_found()
        sys.exit(code)
        
    else:
    
        scaled_mat=motif.getMotif_matrix_scaled()
    
        width=motif.getWidth()
        alphabet=''.join(scaled_mat.index)
        
        pval_mat=np.zeros((width, 100 * width + 1))
    
        for pos in range(width):
            for nuc in alphabet:
                if pos == 0:
                    pval_mat[0, scaled_mat.loc[nuc, pos]] += 1
                else:
                    idxs=np.where(pval_mat[pos - 1, :] > 0)[0]
                    for idx in idxs:
                        pval_mat[pos, scaled_mat.loc[nuc, pos] + idx] += pval_mat[pos-1, idx]
                        
        motif.setMotif_pval_matrix(pval_mat[width-1])
                    
        return motif

def readBGfile(bg_file):
    """
        Read the background file given by the user
        ----
        Parameters:
            bg_file (str) : path to the background file
        ----
        Returns:
            bg_dict (dict) : dictionary that contains the background probabilities
                             defined in the input file
    """
    DNA_ALPHABET=['A', 'C', 'G', 'T']
    
    if bg_file: #we are given the bg file
        
        bg_dict={}
        found_nucs=set()
    
        try:
            bgf=open(bg_file, mode='r') # try to open the file in read only mode
        
            for line in bgf:
                if line[0] in DNA_ALPHABET:
                    nuc, prob=line.split('\t')
                    bg_dict.update({nuc:float(prob)})
                    found_nucs.add(nuc)
                else:
                    # the letter is not part of the DNA alphabet
                    code=he.throw_not_dna_alphabet_error() 
                    sys.exit(code)
                
                if len(found_nucs) == 4:
                    break
                
        except:
            # something went wrong reading the background file
            code=he.throw_not_able_to_read_bg_file()
            sys.exit(code)
        
        finally:
            bgf.close() # close the background file
            
    else: # we are not given of the bg file
        
        bg_dict=get_uniformBG(DNA_ALPHABET)
        
        # assign the uniform distribution to the background
                
    return bg_dict

def get_motif_pwm(motif_file, bgs, pseudo):
    """
        Pipeline to build the motif object
        ----
        Parameters:
            motif_file (str) : path to the motif file
            bgs (str) : path to the background file
            pseudo (float) : value to add to the motif counts
        ----
        Returns:
            motif (Motif) : motif object
    """
    
    if isJaspar_ff(motif_file):
        motif=build_motif_pwm_jaspar(motif_file, bgs, pseudo)
        
    sm=scale_pwm(motif.getMotif_matrix())
    motif.setMotif_matrix_scaled(sm)
        
    motif=comp_pval_mat(motif)
    
    return motif

def get_uniformBG(alphabet):
    """
        Returns a uniform probability distribution for a given alphabet
        ----
        Parameters:
            alphabet (list) : alphabet in list form
        ----
            Returns:
                bg_dict (dict) : dictionary of the uniform probability
    """
    
    if not isinstance(alphabet, list):
        code=he.throw_not_str_error()
        sys.exit(code)
        
    elif isinstance(alphabet, str): # transform the alphabet from a string to
                                    # a list
        alphabet=list(alphabet)
    
    
    alpha_len=len(alphabet)
    bg_dict={}
    up=1/alpha_len
    
    for i in range(alpha_len):
        bg_dict.update({alphabet[i]:up})
            
    return bg_dict
    


                


