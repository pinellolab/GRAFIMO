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
from grafimo.GRAFIMOException import NoDataFrameException, WrongMotifWidthException, WrongMotifIDException, \
                                        WrongMotifNameException, NotValidMotifMatrixException, NotValidBGException, \
                                        NotValidAlphabetException, NotValidFFException, FileReadingException, \
                                        MissingFileException, ValueException
from grafimo.utils import die, DNA_ALPHABET, REV_COMPL, isListEqual, isJaspar_ff, isMEME_ff, almost_equal, \
                            RANGE
from motif_processing import readBGfile, get_uniformBG, apply_pseudocount_jaspar, apply_pseudocount_meme, \
                                compute_log_odds, comp_pval_mat

# define a motif object

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
    
    motif_matrix=None # pd.DataFrame
    motif_scoreMatrix=None # pd.DataFrame
    motif_scoreMatrix_rc=None # pd.DataFrame
    motif_pval_matrix=None # pd.DataFrame
    min_val=-np.inf # double
    max_val=np.inf # double
    scale=-1 # uint
    offset=0 # double
    bg=None # background distribution (dict)
    width=-1
    motifID='' # jaspar ID
    motifName='' # common TF name
    alphabet=None
    isScaled=False
    
    def __init__(self, motif_matrix, width, alphabet, motifID, motifName):
        
        if motif_matrix.empty:
            raise NotValidMotifMatrixException("Attempt to initialize the motif object with an empty motif matrix")
            die(1)
            
        if not isinstance(motif_matrix, pd.DataFrame):
            raise NoDataFrameException("The given value is not a pandas.DatFrame instance")
            die(1)
            
        if not isinstance(width, int) or width < 0:
            raise WrongMotifWidthException("Trying to initialize motif with a not valid width")
            die(1)
            
        if not isinstance(motifID, str) or not motifID:
            raise WrongMotifIDException("Cannot initialize the motif with its ID")
            die(1)

        if not isinstance(motifName, str) or not motifName:
            raise WrongMotifNameException("Cannot initialize the motif with its name")
            die(1)

        if not isinstance(alphabet, list) or not isListEqual(alphabet, DNA_ALPHABET):
            raise NotValidAlphabetException("Cannot initialize a motif object with a wrong alphabet")
            die(1)
        
        self.motif_matrix=motif_matrix
        self.width=width
        self.motifID=motifID
        self.motifName=motifName
        self.alphabet=alphabet

    ### setter methods ###

    def setMotif_matrix(self, motif_matrix):
        
        if motif_matrix.empty:
            raise NotValidMotifMatrixException("Attempt to add to the motif object an empty motif matrix")
            die(1)
            
        if not isinstance(motif_matrix, pd.DataFrame):
            raise NoDataFrameException("The given value is not a pandas.DatFrame instance")
            die(1)

        self.motif_matrix=motif_matrix
            
    def setMotif_scoreMatrix(self, scoreMatrix):
        
        if scoreMatrix.empty:
            raise NotValidMotifMatrixException("Attempt to add to the motif object an empty score matrix")
            die(1)
            
        if not isinstance(scoreMatrix, pd.DataFrame):
            raise NoDataFrameException("The given value is not a pandas.DatFrame instance")
            die(1)

        self.motif_scoreMatrix=scoreMatrix

    def setMotif_scoreMatrix_rc(self, scoreMatrix_rc):

        if scoreMatrix_rc.empty:
            raise NotValidMotifMatrixException("Attempt to add to the motif object an empty score matrix")
            die(1)

        if not isinstance(scoreMatrix_rc, pd.DataFrame):
            raise NoDataFrameException("The given value is not a pandas.DatFrame instance")
            die(1)

        self.motif_scoreMatrix_rc=scoreMatrix_rc
            
    def setMotif_pval_matrix(self, pval_mat):

        if len(pval_mat)==0 or sum(pval_mat[:])<=0: # empty or not valid pvalue matrix
            raise NotValidMotifMatrixException("The p-value matrix is not valid")
        
        self.motif_pval_matrix=pval_mat

    def setMin_val(self, min_val):

        if min_val <= -np.inf:
            raise ValueException(' '.join(["Impossible to assign", min_val, "to Motif.min_val"]))
            die(1)

        self.min_val=min_val

    def setMax_val(self, max_val):

        if max_val >= np.inf:
            raise ValueException(' '.join(["Impossible to assign", max_val, "to Motif.max_val"]))
            die(1)

        self.max_val=max_val

    def setScale(self, scale):

        if not isinstance(scale, int):
            raise ValueException("The scale factor must be an int")
            die(1)

        assert scale > 0

        self.scale=scale

    def setOffset(self, offset):

        self.offset=offset
        
    def setBg(self, bgs):

        if not isinstance(bgs, dict):
            raise NotValidBGException("The background values are not in a dictionary")
            die(1)

        self.bg=bgs
              
    def setWidth(self, width):
        
        if not isinstance(width, int) or width<=0:
            raise WrongMotifWidthException("Trying to initialize motif with a not valid width")
            die(1)

        self.width=width
        
    def setMotifID(self, motifID):
        
        if not isinstance(motifID, str):
            raise WrongMotifIDException(' '.join(["Cannot initialize the motif with the given ID",
                                                    motifID]))
            die(1)
            
        if not motifID:
            raise WrongMotifIDException("Cannot add to the motif with an empty ID")
            die(1)

        self.motifID=motifID
            
    def setMotifName(self, motifName):
        
        if not isinstance(motifName, str):
            raise WrongMotifNameException(' '.join(["Cannot initoalize the motif with the given name",
                                                        motifName]))
            die(1)
            
        if not motifName:
            raise WrongMotifNameException("Cannot add to the motif an empty name")
            die(1)

        self.motifName=motifName
            
    def setAlphabet(self, alphabet):
        
        if not isinstance(alphabet, list):
            raise NotValidAlphabetException("The alphabet given is not in a list")
            die(1)

        if not isListEqual(alphabet, DNA_ALPHABET):
            raise NotValidAlphabetException("The alphabet given is not a valid DNA alphabet")
            die(1)

        self.alphabet=alphabet

    def setIsScaled(self, isScaled):

        if not isinstance(isScaled, bool):
            raise Exception("The isScaled value must be a boolean")
            die(1)

        self.isScaled=isScaled

    ### getter methods ###

    def getMotif_matrix(self):
        
        return self.motif_matrix
    
    def getMotif_scoreMatrix(self):
        
        return self.motif_scoreMatrix

    def getMotif_scoreMatrix_rc(self):

        return self.motif_scoreMatrix_rc
    
    def getMotif_pval_mat(self):
        
        return self.motif_pval_matrix

    def getMin_val(self):

        return self.min_val

    def getMax_val(self):

        return self.max_val

    def getScale(self):

        return self.scale

    def getOffset(self):

        return self.offset
    
    def getBg(self):
        
        return self.bg
    
    def getWidth(self):
        
        return self.width
    
    def getMotifID(self):
        
        return self.motifID
    
    def getMotifName(self):
        
        return self.motifName
    
    def getAlphabet(self):
        
        return self.alphabet

    def getIsScaled(self):

        return self.isScaled
        
    def compute_minValue(self):
        
        motif_matrix=self.getMotif_matrix()
        minValue=motif_matrix.min().sum()
        self.minValue=minValue


# functions to read the motif file and extract from it the score matrix

def build_motif_JASPAR(motif_file, bg_file, pseudocount, no_reverse):
    """
        Build the motif object, given a .jaspar motif file
        ----
        Parameters:
            motif_file (str) : path to the motif file
            bg_file (str) : path to the background file
            pseudocount (float) : value to add to the motif counts
            no_reverse (bool) : boolean value that declares if we want also the
                                reverse complement to be taken into account
        ----
        Returns:
            motif (Motif) : returns a Motif object
    """

    if not motif_file:
        raise MissingFileException("The motif file is missing")
        die(1)

    # check if the input file is in jaspar format
    if not isJaspar_ff(motif_file):
        raise NotValidFFException("The given motif file is not in JASPAR or MEME format")
        die(1)

    assert pseudocount > 0

    # read the motif file
    motif=read_JASPAR_motif(motif_file, bg_file, pseudocount, no_reverse)

    # get the log-odds matrix
    motif_log_odds=compute_log_odds(motif.getMotif_matrix(), motif.getWidth(),
                                        motif.getBg(), motif.getAlphabet())
    motif.setMotif_scoreMatrix(motif_log_odds)

    # scale the log-odds scores
    scaled_scores, min_val, max_val, scale, offset=scale_pwm(motif.getMotif_scoreMatrix(),
                                                                motif.getAlphabet(),
                                                                motif.getWidth())
    motif.setMotif_scoreMatrix(scaled_scores)
    motif.setIsScaled(True)
    motif.setScale(scale)
    motif.setMin_val(min_val)
    motif.setMax_val(max_val)
    motif.setOffset(offset)

    # compute the p-value matrix
    pval_mat=comp_pval_mat(motif)
    motif.setMotif_pval_matrix(pval_mat)
            
    return motif

def read_JASPAR_motif(motif_file, bg_file, pseudocount, no_reverse):
    """
        Read the motif file in JASPAR format and build a motif
        object from it
        ----
        Params:
            motif_file (str) : path to the motif file (in JASPAR format)
            bg_file (str) : path to the background file
            no_reverse (bool) : boolean value that declares if we want also the
                                reverse complement to be taken into account
        ----
        Returns:
            motif (Motif) : motif object built from the data in motif_file
    """

    nucs=[]
    counts=[]

    try:
        mf=open(motif_file, mode='r')  # open the file in read only mode

        header=str(mf.readline()[1:])  # read the header
        motifID, motifName = header.split('\t')[0:2]  # get the jaspar ID and the common TF name
        motifName=motifName[:-1]  # avoid '\n'

        for line in mf:
            line=line.strip()
            nuc=line.strip()[:1]
            count=list(map(float, line.strip()[1:].split()[1:][:-1]))

            nucs.append(nuc)
            counts.append(count)

    except:  # something went wrong
        msg=' '.join(["Unable to read file", motif_file])
        raise FileReadingException(msg)
        die(1)

    else:

        motif_counts=pd.DataFrame(data=counts, index=nucs)  # raw counts
        motif_width=int(len(motif_counts.columns))
        alphabet=sorted(nucs)

        if bg_file:
            bgs=readBGfile(bg_file)
        else:
            bgs=get_uniformBG(alphabet)

        bgs=pseudo_bg(bgs, no_reverse)

        motif_probs=(motif_counts / motif_counts.sum(0))  # get probabilities
        motif_probs=norm_motif(motif_probs, motif_width, alphabet)
        motif_probs=apply_pseudocount_jaspar(motif_counts, motif_probs, pseudocount, bgs,
                                              motif_width, alphabet)

        motif=Motif(motif_probs, motif_width, alphabet, motifID, motifName)
        motif.setBg(bgs)

        return motif

    finally:
        mf.close() # close the motif file anyway


def build_motif_MEME(motif_file, bg_file, pseudocount, no_reverse):
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

    if not motif_file:
        raise MissingFileException("The motif file is missing")
        die(1)

    # check if the input is in MEME format
    if not isMEME_ff(motif_file):
        # if in other format we should not be here
        raise Exception("Error: the motif file given is not in MEME format") 
        sys.exit(1)

    # read the motif file
    motif=read_MEME_motif(motif_file, bg_file, pseudocount, no_reverse)

    # get the log-odds
    motif_log_odds=compute_log_odds(motif.getMotif_matrix(), motif.getWidth(),
                                        motif.getBg(), motif.getAlphabet())
    motif.setMotif_scoreMatrix(motif_log_odds)

    # scale the log-odds scores
    scaled_scores, min_val, max_val, scale, offset=scale_pwm(motif.getMotif_scoreMatrix(),
                                                                motif.getAlphabet(),
                                                                motif.getWidth())
    motif.setMotif_scoreMatrix(scaled_scores)
    motif.setIsScaled(True)
    motif.setScale(scale)
    motif.setMin_val(min_val)
    motif.setMax_val(max_val)
    motif.setOffset(offset)

    # compute the p-value matrix
    pval_mat=comp_pval_mat(motif)
    motif.setMotif_pval_matrix(pval_mat)

    return motif


def read_MEME_motif(motif_file, bg_file, pseudocount, no_reverse):
    """
        Read the motif file in MEME format and build a motif
        object from it
        ----
        Params:
            motif_file (str) : path to the motif file
            bg_file (str) : path to the background file
            pseudocount (np.double) : pseudocount to apply to the motif
            no_reverse (bool) : boolean value that declares if we wnat also the
                                reverse complement to be taken into account
        ----
        Returns:
            motif (Motif) returns a Motif object
    """

    try:
        mf = open(motif_file, 'r')  # open the file in read only mode

        infostart = False  # flag to keep track were the infos about the motif begin
        datastart = False  # flag to keep track were the motif data begin

        a = []
        c = []
        g = []
        t = []

        pos_read=0

        for line in mf:
            if line[0:8]=='ALPHABET':
                alphabet = sorted(list(set(line[10:-1])))
                assert isListEqual(alphabet, DNA_ALPHABET)

            if line[0:5]=='MOTIF':
                motifID, motifName = line.split(' ')[1:3]
                motifName=motifName[:-1]  # remove \n char
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
                infostart=False  # informations end here

                # allocate space for the motif probability matrix
                motif_probs=pd.DataFrame(index=alphabet, columns=range(motif_width),
                                          data=np.double(0))

                datastart = True  # at next step begin data
                continue

            if datastart and pos_read < motif_width:
                freqs=line.split(' ')
                a.append(np.double(freqs[1]))
                c.append(np.double(freqs[3]))
                g.append(np.double(freqs[5]))
                t.append(np.double(freqs[7][:-1]))  # remove '\n'
                pos_read+=1

    except:  # something went wrong
        msg = ' '.join(["Unable to read file", motif_file])
        raise FileReadingException(msg)
        die(1)

    else:
        motif_probs.loc['A']=a
        motif_probs.loc['C']=c
        motif_probs.loc['G']=g
        motif_probs.loc['T']=t

        if bg_file:
            bgs=readBGfile(bg_file)
        else:
            bgs=get_uniformBG(alphabet)

        bgs=pseudo_bg(bgs, no_reverse)

        motif_probs=norm_motif(motif_probs, motif_width, alphabet)
        motif_probs=apply_pseudocount_meme(motif_probs, pseudocount, site_counts, motif_width,
                                            bgs, alphabet)

        motif=Motif(motif_probs, motif_width, alphabet, motifID, motifName)
        motif.setBg(bgs)

        return motif

    finally:
        mf.close() # close the file anyway


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
        raise NoDataFrameException("The given motif matrix must be an instance of pandas.DataFrame")
        die(1)

    if motif_matrix.empty:
        raise NotValidMotifMatrixException("The given motif matrix is empty")
        die(1)

    if not isinstance(alphabet, list):
        raise NotValidAlphabetException("The alphabet given is not in a list")
        die(1)

    if not isListEqual(alphabet, DNA_ALPHABET):
        raise NotValidAlphabetException("The alphabet given is not a valid DNA alphabet")
        die(1)

    assert motif_width > 0

    min_val=min(motif_matrix.min())
    max_val=max(motif_matrix.max())
    motif_matrix_sc=pd.DataFrame(index=list(motif_matrix.index), columns=list(motif_matrix.columns),
                                        data=0)

    lower=min_val
    upper = max_val
        
    if lower == upper: # all values are equal
        lower=np.double(upper-1)
        
    lower=np.floor(lower)
    offset=np.round(np.floor(lower))
    scale_factor=np.floor(RANGE/(upper-lower))

    # values will be in [0, 1000]
    for nuc in alphabet:
        for j in range(motif_width):
            scaled_score=np.round((motif_matrix.loc[nuc, j]-(offset))*scale_factor)
            motif_matrix_sc.loc[nuc, j]=scaled_score
    
    # make sure the values are integers
    motif_matrix_sc[:]=motif_matrix_sc[:].astype(int)

    # now they are scaled
    min_val = min(motif_matrix_sc.min())
    max_val = max(motif_matrix_sc.max())
     
    return motif_matrix_sc, min_val, max_val, int(scale_factor), offset


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

    if not motif_file:
        raise MissingFileException("The motif file is missing")
        die(1)

    if isJaspar_ff(motif_file):
        motif=build_motif_JASPAR(motif_file, bgs, pseudo, no_reverse)
        
    elif isMEME_ff(motif_file):
        motif=build_motif_MEME(motif_file, bgs, pseudo, no_reverse)

    else:
        msg=' '.join["Do not know what to do with file", motif_file]
        raise NotValidFFException(msg)
        die(1)

    return motif

def pseudo_bg(bgs, no_reverse):
    """ 
        Apply a pseudo count to the background frequencies
        ----
        Parameters:
            bgs (dict) : dictionary of the background probabilities
        ----
        Returns:
            bgs_proc (dict) : normalized (and averaged) background frequencies
    """
        
    if not isinstance(bgs, dict):
        raise NotValidBGException("Failed attempt to add pseudocount to the background")
        die(1)

    if not isinstance(no_reverse, bool):
        raise ValueException(' '.join(["Boolean value required, got", str(type(no_reverse))]))
        die(1)
        
    if not no_reverse:
        bgs_avg=average_bg_with_rc(bgs)
    else:
        bgs_avg=bgs
        
    bgs_proc=norm_bg(bgs_avg)

    return bgs_proc
        
def average_bg_with_rc(bgs):
    """
        Average the background frequencies with the reverse complement
        ----
        Parameters:
            bgs (dict) : dictionary of the background probabilities
        ----
        Returns:
            bgs_avg (dict) : averaged background probabilities
    """

    bgs_avg={}
    
    for nuc in bgs.keys():
        rc=REV_COMPL[nuc]

        if REV_COMPL[rc]==nuc and ord(nuc) < ord(rc):
            avg_freq=np.double((bgs[nuc]+bgs[rc])/np.double(2))
            bgs_avg.update({nuc: avg_freq})
            bgs_avg.update({rc:avg_freq})
            
    return bgs_avg

def norm_bg(bgs):
    """
        Normalize the bakground frequencies
        ----
        Parameters:
            Parameters:
            bgs (dict) : dictionary of the background probabilities
        ----
        Returns:
            bgs_norm (dict) : normalized background probabilities
    """
    
    PSEUDO=np.double(0.0000005) # pseudocount

    alphabet=sorted(list(bgs.keys()))
    tot=np.double(len(alphabet)*PSEUDO)
    bgs_norm={}
    
    for nuc in bgs.keys():
        tot+=np.double(bgs[nuc])

    assert tot > 0
        
    for nuc in bgs.keys():
        prob=np.double((bgs[nuc]+PSEUDO)/tot)
        bgs_norm.update({nuc : prob})
        
    tot=np.double(0)
    for nuc in bgs.keys():
        tot+=bgs[nuc]
    
    assert tot!=0
    
#    if almost_equal(1, tot, 0):
#        print("is almost equal")
#        return bgs_norm
#    else:
#        for nuc in bgs.keys():
#            bgs_norm[nuc]=np.double(bgs_norm[nuc]/tot)
            
    return bgs_norm

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
                      # probability and 1
                      
    for j in range(motif_width):
        tot=np.double(0)

        for nuc in alphabet:
            tot+=motif_probs.loc[nuc, j]

        assert tot != 0
            
        if not almost_equal(1, tot, tolerance):
            for nuc in alphabet:
                motif_probs.loc[nuc, j]=np.double(motif_probs.loc[nuc, j]/tot)
                
    return motif_probs

