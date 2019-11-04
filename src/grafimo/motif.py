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
    bg={} # background distribution
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
        
    def setBg(self, bgs):
        if not isinstance(bgs, dict):
            raise Exception("Error: unable to set the background distribution")
            sys.exit(1)
        
        else:
            self.bg=bgs
              
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
    
    def getBg(self):
        
        return self.bg
    
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
    
def isMEME_ff(motif_file):
    """
        Check if the given file is in .meme format
        ----
        Parameters:
            motif_file (str) : path to the motif file
        ----
        Returns:
            (bool)
    """
    
    if motif_file and isinstance(motif_file, str):
        ff=motif_file.split('.')[-1]
        
        if ff=='meme':
            return True
        else:
            return False
        
    else:
        return False # the motif file was not given or the path is empty
            
def build_motif_pwm_jaspar(motif_file, bg_file, pseudocount, no_reverse):
    """
        Build the motif object, given a .jaspar motif file
        ----
        Parameters:
            motif_file (str) : path to the motif file
            bg_file (str) : path to the background file
            pseudocount (float) : value to add to the motif counts
            no_reverse (bool) : boolean value that declares if we wnat also the
                                reverse complement to be taken into account
        ----
        Returns:
            motif (Motif) : returns a Motif object
    """
    
    # check if the input file is in jaspar format
    if not isJaspar_ff(motif_file):
        raise Exception("Error: the given mnotif file is not in JASPAR format")
        sys.exit(1)
    
    # build the motif pwm
    nucs = []
    fqs = []
    
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
        msg="Error: uanble to read motif "+motif_file
        raise Exception(msg)
        sys.exit(1)
            
    else:    
            
        motif_count = pd.DataFrame(data=fqs, index=nucs) # raw counts
        motif_width=len(motif_count.columns)
        bgs=readBGfile(bg_file)
        alphabet=sorted(nucs)
        
        bgs=pseudo_bg(bgs, no_reverse)
    
        motif_prob=(motif_count/motif_count.sum(0)) # get probabilities
        motif_prob=norm_motif(motif_prob, motif_width, alphabet)
   
        for j in range(motif_width): 
            site_counts=sum(motif_count.loc[:,j])
            total_counts=site_counts+pseudocount
            for nuc in alphabet:
                bg=bgs[nuc]
                count=np.double((motif_prob.loc[nuc, j]*site_counts)+
                                    (pseudocount*bg))
                motif_prob.loc[nuc, j]=np.double(count/total_counts)
                
        # allocate the log-odds matrix
        motif_log_odds=pd.DataFrame(index=nucs, columns=range(motif_width), 
                                        data=np.double(0))
            
        tot_bg=np.double(0)
        tot_fg=np.double(0)
        
        for nuc in alphabet:
            bg=bgs[nuc]
            tot_bg+=bg
            for j in range(motif_width):
                prob=motif_prob.loc[nuc, j]
                tot_fg+=prob
                logodds=np.double(np.log2(prob/bg))
                motif_log_odds.loc[nuc, j]=logodds
                
        # check if both the background and the foreground are next to 1
        epsilon=0.001
        assert tot_bg-1<epsilon
        assert tot_fg-motif_width<epsilon
            
        motif=Motif(motif_log_odds, motif_width, motifID=motifID, 
                        motifName=motifName) # create the motif
        alphabet=''.join(alphabet)
        motif.setAlphabet(alphabet)
        motif.setBg(bgs)
            
        return motif
    
    finally:
        mf.close() # close the motif file anyway

def build_motif_pwm_MEME(motif_file, bg_file, pseudocount, no_reverse):
    """
        Build the motif object, given a .jaspar motif file
        ----
        Parameters:
            motif_file (str) : path to the motif file
            bg_file (str) : path to the background file
            pseudocount (float) : value to add to the motif counts
            no_reverse (bool) : boolean value that declares if we wnat also the
                                reverse complement to be taken into account
        ----
        Returns:
            motif (Motif) : returns a Motif object
    """
    
    if not isMEME_ff(motif_file):
        # if in other format we should not be here
        raise Exception("Error: the motif file given is not in MEME format") 
        sys.exit(1)
        
    try:
        mf=open(motif_file, 'r') # open the file in read only mode
        
        infostart=False # flag to keep track were the infos about the motif begin
        datastart=False # flag to keep track were the motif data begin
        
        A=[]
        C=[]
        G=[]
        T=[]
        
        pos_read=0
        
        for line in mf:
            if line[0:8]=='ALPHABET':
                alphabet=sorted(list(set(line[10:-1])))
                
            if line[0:5]=='MOTIF':
                motifID, motifName=line.split[1:3]
                motifName=motifName[:-1] # remove \n char
                # the informations about motif start here
                infostart=True
                continue
            
            if infostart:
                infos=line[26:]
                infosplit=infos.split(' ')
                alphalen=int(infosplit[2])
                
                assert alphalen==len(alphabet)
                
                motif_width=int(infosplit[4])
                site_counts=int(infosplit[6])
                infostart=False # informations end here
                
                # allocate space for the motif probability matrix
                motif_prob=pd.DataFrame(index=alphabet, 
                                            columns=range(motif_width),
                                            data=np.double(0))
                
                datastart=True # at next step begin data
                continue
            
            if datastart and pos_read<motif_width:
                freqs=line.split(' ')
                A.append(np.double(freqs[1]))
                C.append(np.double(freqs[3]))
                G.append(np.double(freqs[5]))
                T.append(np.double(freqs[7][:-1])) # remove \n char
                pos_read+=1
                
    except: # something went wrong
        msg="Error: unable to read motif "+motif_file
        raise Exception(msg)
        sys.exit(1)
        
    else:
        motif_prob.loc['A']=A
        motif_prob.loc['C']=C
        motif_prob.loc['G']=G
        motif_prob.loc['T']=T
        
        bgs=readBGfile(bg_file)
        
        bgs=pseudo_bg(bgs, no_reverse) 
        
        motif_prob=norm_motif(motif_prob, motif_width, alphabet)
        
        # build the scoring matrix
        total_counts=site_counts+pseudocount
        
        for j in range(motif_width):
            for nuc in alphabet:
                bg=bgs[nuc]
                count=np.double((motif_prob.loc[nuc, j]*site_counts)+
                                    (pseudocount*bg))
                motif_prob.loc[nuc, j]=np.double(count/total_counts)
                
        motif_log_odds=pd.DataFrame(index=alphabet, columns=range(motif_width),
                                        data=np.double(0))
        tot_bg=np.double(0)
        tot_fg=np.double(0)
        
        for nuc in alphabet:
            bg=bgs[nuc]
            tot_bg+=bg
            for j in range(motif_width):
                prob=motif_prob.loc[nuc, j]
                tot_fg+=prob
                logodds=np.double(np.log2(prob/bg))
                motif_log_odds.loc[nuc, j]=logodds
                
        # check if both the background and the foreground are next to 1
        epsilon=0.001
        assert tot_bg-1<epsilon
        assert tot_fg-motif_width<epsilon
        
        motif=Motif(motif_log_odds, motif_width, motifID=motifID, 
                        motifName=motifName)
        alphabet=''.join(alphabet)
        motif.setAlphabet(alphabet)
        motif.setBg(bgs)
        
        return motif
    
    finally:
        mf.close() # close the motif file stream anyway
        
def scale_pwm(motif_matrix, alphabet, motif_width):
    """
        Computes the scaled matrix, given a motif matrix
        ----
        Parameters:
            motif_matrix (str) : pat to the motif file
            alphabet (str) : alphabet of the motif
            motif_width (int) : width of the motif
        ----
        Returns:
            motif_matrix (pd.DataFrame) : scaled motif matrix
    """
    
    if not isinstance(motif_matrix, pd.DataFrame):
        raise Exception("Error: unable to read the motif scoring matrix")
        sys.exit(1)
        
    else:
    
        RANGE_VALUE=1000 # values will be in [0, 1000]
        min_val=min(motif_matrix.min())
        max_val=max(motif_matrix.max())
        motif_matrix_sc=pd.DataFrame(index=list(motif_matrix.index),
                                         columns=list(motif_matrix.columns),
                                         data=0)
        
        if min_val==max_val: # all values are equal
            min_val=max_val-1
        
        offset=np.floor(min_val)
        scale_factor=np.floor(RANGE_VALUE/(max_val-min_val))
        
        for nuc in alphabet:
            for j in range(motif_width):
                scaled_score=np.round((motif_matrix.loc[nuc, j]-(1*offset))*scale_factor)
                motif_matrix_sc.loc[nuc, j]=scaled_score
    
        # make sure the values are integers
        motif_matrix_sc[:]=motif_matrix_sc[:].astype(int)
     
        return motif_matrix_sc

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
        msg="Error: the given object is not an instance of Motif. Unable to compute the p-value matrix"
        raise Exception(msg)
        sys.exit(1)
        
    if not motif.has_scaled_matrix():
        msg="Error: the scaled scoring matrix is missing. Unable to compute the p-value matrix"
        raise Exception(msg)
        sys.exit(1)
        
    RANGE=1000
    scaled_mat=motif.getMotif_matrix_scaled()
    
    motif_width=motif.getWidth()
    alphabet=motif.getAlphabet()
    bgs=motif.getBg()
        
    pval_mat=np.zeros((motif_width, RANGE*motif_width + 1))
    
    for pos in range(motif_width):
        for nuc in alphabet:
            if pos==0:
                pval_mat[0, scaled_mat.loc[nuc, pos]] += 1
            else:
                idxs=np.where(pval_mat[pos - 1, :]>0)[0]
                bg=bgs[nuc]
                for idx in idxs:
                    pval_mat[pos, scaled_mat.loc[nuc, pos] + idx]+=pval_mat[pos-1, idx]*bg
                        
    motif.setMotif_pval_matrix(pval_mat[motif_width-1])
                    
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
            raise Exception("Error: unable to read the background file")
            sys.exit(1)
        
        finally:
            bgf.close() # close the background file
            
    else: # we are not given of the bg file
        
        bg_dict=get_uniformBG(DNA_ALPHABET)
        
        # assign the uniform distribution to the background
                
    return bg_dict

def get_motif_pwm(motif_file, bgs, pseudo, no_reverse):
    """
        Pipeline to build the motif object
        ----
        Parameters:
            motif_file (str) : path to the motif file
            bgs (str) : path to the background file
            pseudo (float) : value to add to the motif counts
            no_reverse (bool) : flag that says if has to be taken into account 
                                also the reverse complement
        ----
        Returns:
            motif (Motif) : motif object
    """
    
    if isJaspar_ff(motif_file):
        motif=build_motif_pwm_jaspar(motif_file, bgs, pseudo, no_reverse)
        
    elif isMEME_ff(motif_file):
        motif=build_motif_pwm_MEME(motif_file, bgs, pseudo, no_reverse)
        
    sm=scale_pwm(motif.getMotif_matrix(), motif.getAlphabet(), 
                    motif.getWidth())

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

def pseudo_bg(bgs, no_reverse):
    """ 
        Apply the pseudo count to the background frequencies
        ----
        Parameters:
            bgs (dict) : dictionary of the background probabilities
        ----
        Returns:
            bgs (dict) : normalized (and normalized) background frequencies
    """
        
    if not isinstance(bgs, dict):
        raise Exception('Error: unable to normalize background frequencies')
        sys.exit(1)
        
    if not no_reverse:
        bgs=average_bg_with_rc(bgs)
        
    bgs=norm_bg(bgs)
    
    return bgs
        
def average_bg_with_rc(bgs):
    """
        Average the background frequencies with the reverse complement
        ----
        Parameters:
            bgs (dict) : dictionary of the background probabilities
        ----
        Returns:
            bgs (dict) : averaged background probabilities
    """
    
    RC={'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    
    for nuc in bgs.keys():
        rc=RC[nuc]
        if RC[rc]==nuc and ord(nuc)<ord(rc):
            avg_freq=np.double((bgs[nuc]+bgs[rc])/2.0)
            bgs[nuc]=avg_freq
            bgs[rc]=avg_freq
            
    return bgs

def norm_bg(bgs):
    """
        Normalize the bakground frequencies
        ----
        Parameters:
            Parameters:
            bgs (dict) : dictionary of the background probabilities
        ----
        Returns:
            bgs (dict) : normalized background probabilities
    """
    
    pseudo=np.double(0.0000005)
    alphabet=sorted(list(bgs.keys()))
    tot=np.double(len(alphabet)*pseudo)
    
    for nuc in bgs.keys():
        tot+=bgs[nuc]
        
    for nuc in bgs.keys():
        bgs[nuc]=np.double((bgs[nuc]+pseudo)/tot)
        
    tot=np.double(0)
    for nuc in bgs.keys():
        tot+=bgs[nuc]
    
    assert tot!=0
    
    if not almost_equal(1, tot, 0):
        for nuc in bgs.keys():
            bgs[nuc]=np.double(bgs[nuc]/tot)
            
    return bgs

def norm_motif(motif_probs, motif_width, alphabet):
    """
        Normalize the motif probabilities
        ----
        Parameters:
            motif_probs (pd.DataFarme) : probability matrix
            motif_width (int) : width of the motif
            alphabet (list) : alphabet of the motif
        ----
        Returns:
            motif_probs (pd.DataFrame) : normalized probability matrix
    """
    
    tolerance=0.00001 # allowed tolerance in the difference between the position
                      # probaility and 1
                      
    for j in range(motif_width):
        tot=np.double(0)
        for nuc in alphabet:
            tot+=motif_probs.loc[nuc, j]
            
        if not almost_equal(1, tot, tolerance):
            for nuc in alphabet:
                motif_probs.loc[nuc, j]=np.double(motif_probs.loc[nuc, j]/tot)
                
    return motif_probs
        
        
def almost_equal(value1, value2, slope):
    
    if (value1-slope)>value2 or (value1+slope)<value2:
        return False
    else:
        return True

