"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

The file contains the definition of the Motif class, that represent
the motif that will be searched on the genome.

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


#################################################
# begin Motif definition
#################################################
class Motif(object):
    """
        Class to represent a motif.
        It contains:
            - the motif count matrix
            - the motif scoring matrix (scaled)
            - the pvalue matrix
            - parameters used to scale the matrix
            - the background distribution of the motif
            - the width of the motif
            - the minimum value in the motif matrix (both forward and reverse)
            - the maximum value in the motif matrix (both forward and reverse)
            - the transcription factor's name (both jaspar ID and conventional name)
            - the alphabet on which the matrix is built
    """
    
    motif_matrix = None # pd.DataFrame
    motif_scoreMatrix = None # np.ndarray
    motif_pval_matrix = None # np.array
    min_val = -np.inf # np.double
    max_val = np.inf # np.double
    scale = -1 # int
    offset = 0 # np.double
    bg = None # background distribution (dict)
    width = -1 # int
    motifID = '' # jaspar ID
    motifName = '' # common TF name
    alphabet = None
    isScaled = False
    
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
        
        self.motif_matrix = motif_matrix
        self.width = width
        self.motifID = motifID
        self.motifName = motifName
        self.alphabet = alphabet

    ### setter methods ###

    def setMotif_matrix(self, motif_matrix):
        
        if motif_matrix.empty:
            raise NotValidMotifMatrixException("Attempt to add to the motif object an empty motif matrix")
            die(1)
            
        if not isinstance(motif_matrix, pd.DataFrame):
            raise NoDataFrameException("The given value is not a pandas.DatFrame instance")
            die(1)

        self.motif_matrix = motif_matrix
            
    def setMotif_scoreMatrix(self, scoreMatrix):
            
        if not isinstance(scoreMatrix, np.ndarray) and \
            not isinstance(scoreMatrix, pd.DataFrame):
            raise ValueError("The given data structure is not an instance of numpy.ndarray or pandas.DataFrame")
            die(1)

        if isinstance(scoreMatrix, pd.DataFrame):
            if scoreMatrix.empty:
                raise NotValidMotifMatrixException("Attempt to add to the motif object an empty score matrix")
                die(1)

        if isinstance(scoreMatrix, np.ndarray):
            if scoreMatrix.size == 0:
                raise NotValidMotifMatrixException("Attempt to add to the motif object an empty score matrix")
                die(1)

        self.motif_scoreMatrix = scoreMatrix
            
    def setMotif_pval_matrix(self, pval_mat):

        if len(pval_mat) == 0 or sum(pval_mat[:]) <= 0: # empty or not valid p-value matrix
            raise NotValidMotifMatrixException("The p-value matrix is not valid")
        
        self.motif_pval_matrix = pval_mat

    def setMin_val(self, min_val):

        if min_val <= -np.inf:
            raise ValueException(' '.join(["Impossible to assign", min_val, "to Motif.min_val"]))
            die(1)

        self.min_val = min_val

    def setMax_val(self, max_val):

        if max_val >= np.inf:
            raise ValueException(' '.join(["Impossible to assign", max_val, "to Motif.max_val"]))
            die(1)

        self.max_val = max_val

    def setScale(self, scale):

        if not isinstance(scale, int):
            raise ValueException("The scale factor must be an int")
            die(1)

        assert scale > 0

        self.scale = scale

    def setOffset(self, offset):

        self.offset = offset
        
    def setBg(self, bgs):

        if not isinstance(bgs, dict):
            raise NotValidBGException("The background values are not in a dictionary")
            die(1)

        self.bg = bgs
              
    def setWidth(self, width):
        
        if not isinstance(width, int) or width <= 0:
            raise WrongMotifWidthException("Trying to initialize motif with a not valid width")
            die(1)

        self.width = width
        
    def setMotifID(self, motifID):
        
        if not isinstance(motifID, str):
            raise WrongMotifIDException(' '.join(["Cannot initialize the motif with the given ID",
                                                    motifID]))
            die(1)
            
        if not motifID:
            raise WrongMotifIDException("Cannot add to the motif with an empty ID")
            die(1)

        self.motifID = motifID
            
    def setMotifName(self, motifName):
        
        if not isinstance(motifName, str):
            raise WrongMotifNameException(' '.join(["Cannot initoalize the motif with the given name",
                                                        motifName]))
            die(1)
            
        if not motifName:
            raise WrongMotifNameException("Cannot add to the motif an empty name")
            die(1)

        self.motifName = motifName
            
    def setAlphabet(self, alphabet):
        
        if not isinstance(alphabet, list):
            raise NotValidAlphabetException("The alphabet given is not in a list")
            die(1)

        if not isListEqual(alphabet, DNA_ALPHABET):
            raise NotValidAlphabetException("The alphabet given is not a valid DNA alphabet")
            die(1)

        self.alphabet = alphabet

    def setIsScaled(self, isScaled):

        if not isinstance(isScaled, bool):
            raise Exception("The isScaled value must be a boolean")
            die(1)

        self.isScaled = isScaled

    ### getter methods ###

    def getMotif_matrix(self):
        
        return self.motif_matrix
    
    def getMotif_scoreMatrix(self):
        
        return self.motif_scoreMatrix
    
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

#################################################
# end of Motif definition
#################################################


def build_motif_JASPAR(motif_file, bg_file, pseudocount, no_reverse):
    """
        Build a the Motif object starting from the raw counts
        data stored in a given JASPAR file.

        The raw counts are processed and the resulting values 
        are used to define the scoring matrix of the motif
        ----
        Parameters:
            motif_file (str) : path to the motif file
            bg_file (str) : path to the background file
            pseudocount (float) : value to add to the motif counts (to avoid 
                                    division by 0)
            no_reverse (bool) : flag parameter to consider or not the reverse 
                                complement building the Motif object
        ----
        Returns:
            motif (Motif) : returns the corresponding Motif object
    """

    # This shouldn't happen
    if not motif_file:
        raise MissingFileException("The motif file is missing")
        die(1)

    # check if the input file is in jaspar format
    if not isJaspar_ff(motif_file):
        raise NotValidFFException("The given motif file is not in JASPAR or MEME format")
        die(1)

    assert pseudocount > 0

    # read the motif file
    motif = read_JASPAR_motif(motif_file, bg_file, pseudocount, no_reverse)

    # get the log-odds matrix
    motif_log_odds = compute_log_odds(motif.getMotif_matrix(), motif.getWidth(),
                                        motif.getBg(), motif.getAlphabet())
    motif.setMotif_scoreMatrix(motif_log_odds)

    # scale the log-odds scores
    scaled_scores, min_val, max_val, scale, offset = scale_pwm(motif.getMotif_scoreMatrix(),
                                                                motif.getAlphabet(),
                                                                motif.getWidth())

    motif.setMotif_scoreMatrix(scaled_scores)
    motif.setIsScaled(True)
    motif.setScale(scale)
    motif.setMin_val(min_val)
    motif.setMax_val(max_val)
    motif.setOffset(offset)

    # compute the p-value matrix
    pval_mat = comp_pval_mat(motif)
    motif.setMotif_pval_matrix(pval_mat)

    motif.setMotif_scoreMatrix(scaled_scores.values) # store the np.ndarray of the score matrix
            
    return motif


def read_JASPAR_motif(motif_file, bg_file, pseudocount, no_reverse):
    """
        Read the data contained in a JASPAR file and build a motif
        object from it
        ----
        Params:
            motif_file (str) : path to the motif file (in JASPAR format)
            bg_file (str) : path to the background file
            no_reverse (bool) : flag parameter to consider or not the reverse 
                                complement building the Motif object
        ----
        Returns:
            motif (Motif) : Motif object built from the data in motif_file
    """

    nucs = []
    counts = []

    try:
        mf = open(motif_file, mode='r')  # open the file in read only mode

        header = str(mf.readline()[1:])  # read the header
        motifID, motifName = header.split('\t')[0:2]  # get the jaspar ID and the common TF name
        motifName = motifName[:-1]  # avoid '\n'

        for line in mf:
            line = line.strip()
            nuc = line.strip()[:1]
            count = list(map(float, line.strip()[1:].split()[1:][:-1]))

            nucs.append(nuc)
            counts.append(count)

    except:  # something went wrong
        msg = ' '.join(["Unable to read file", motif_file])
        raise FileReadingException(msg)
        die(1)

    else:

        motif_counts = pd.DataFrame(data=counts, index=nucs)  # raw counts
        motif_width = int(len(motif_counts.columns))
        alphabet = sorted(nucs)

        if bg_file:
            bgs = readBGfile(bg_file)
        else:
            bgs = get_uniformBG(alphabet)

        bgs = pseudo_bg(bgs, no_reverse)

        motif_probs = (motif_counts / motif_counts.sum(0))  # get probabilities
        motif_probs = norm_motif(motif_probs, motif_width, alphabet)
        motif_probs = apply_pseudocount_jaspar(motif_counts, motif_probs, pseudocount, bgs,
                                                motif_width, alphabet)

        motif = Motif(motif_probs, motif_width, alphabet, motifID, motifName)
        motif.setBg(bgs)

        return motif

    finally:
        mf.close() # close the motif file anyway


def build_motif_MEME(motif_file, bg_file, pseudocount, no_reverse):
    """
        Build a the Motif object starting from the probabilities
        stored in a given MEME file.

        The probabilities are processed and the resulting values 
        are used to define the scoring matrix of the motif.
        ----
        Parameters:
            motif_file (str) : path to the motif file
            bg_file (str) : path to the background file
            pseudocount (float) : value to add to the motif counts
            no_reverse (bool) : flag parameter to consider or not the reverse 
                                complement building the Motif object
        ----
        Returns:
            motif (Motif) : Motif object built from data contained in
                            motif_file
    """

    # This shouldn't happen
    if not motif_file:
        raise MissingFileException("The motif file is missing")
        die(1)

    # check if the input is in MEME format
    if not isMEME_ff(motif_file):
        # if in other format we should not be here
        raise NotValidFFException("Error: the motif file given is not in MEME format") 
        sys.exit(1)

    # read the motif file
    motif_lst = read_MEME_motif(motif_file, bg_file, pseudocount, no_reverse)
    
    # list of the fully processed motifs
    complete_motifs = []

    # process each found motif
    for m in motif_lst:
        
        # get the log-odds
        motif_log_odds = compute_log_odds(m.getMotif_matrix(), m.getWidth(),
                                                m.getBg(), m.getAlphabet())
        m.setMotif_scoreMatrix(motif_log_odds)

        # scale the log-odds scores
        scaled_scores, min_val, max_val, scale, offset = scale_pwm(m.getMotif_scoreMatrix(),
                                                                    m.getAlphabet(),
                                                                    m.getWidth())
        m.setMotif_scoreMatrix(scaled_scores)
        m.setIsScaled(True)
        m.setScale(scale)
        m.setMin_val(min_val)
        m.setMax_val(max_val)
        m.setOffset(offset)

        # compute the p-value matrix
        pval_mat = comp_pval_mat(m)
        m.setMotif_pval_matrix(pval_mat)

        m.setMotif_scoreMatrix(scaled_scores.values) 

        complete_motifs.append(m)

    return complete_motifs


def read_MEME_motif(motif_file, bg_file, pseudocount, no_reverse):
    """
        Read the motif file in MEME format and build a motif
        object from it.
        Note that a MEME file can contain a variable number of
        motifs
        ----
        Params:
            motif_file (str) : path to the motif file
            bg_file (str) : path to the background file
            pseudocount (np.double) : pseudocount to apply to the motif
            no_reverse (bool) : flag parameter to consider or not the reverse 
                                complement building the Motif object
        ----
        Returns:
            motif (Motif) : returns a Motif object
    """

    try:
        mf = open(motif_file, 'r')  # open the file in read only mode

        infostart = False  # flag to keep track were the infos about the motif begin
        datastart = False  # flag to keep track were the motif data begin
        motifs_found = 0 # number of motifs found in the MEME file 
        
        motifID_lst = [] # list of the found motif IDs
        motifName_lst = [] # list of the found motif names
        motif_width_lst = [] # list of the found motif widths
        site_counts_lst = [] # list of the found motif site counts
        alphalen_lst = [] # list of the found motif alphabet lengths
        motif_probs_lst = [] # list of the found motif probability matrices
        a_lst = [] # list of the found As probabilities for each motif
        c_lst = [] # list of the found Cs probabilities for each motif
        g_lst = [] # list of the found Gs probabilities for each motif
        t_lst = [] # list of the found Ts probabilities for each motif

        motif_width = None
        pos_read = 0

        for line in mf:
            if line[0:8] == 'ALPHABET':
                alphabet = sorted(list(set(line[10:-1])))
                assert isListEqual(alphabet, DNA_ALPHABET)

            if line[0:5] == 'MOTIF':
                motifID, motifName = line.split()[1:3]

                motifID_lst.append(motifID)
                motifName_lst.append(motifName)

                # the informations about motif start here
                infostart = True
                continue

            if infostart and len(line.strip()) != 0:
                infos = line[26:]
                infosplit = infos.split()
                alphalen = int(infosplit[1])
                alphalen_lst.append(alphalen)

                assert alphalen == len(alphabet)

                motif_width = int(infosplit[3])
                site_counts = int(infosplit[5])
                infostart = False  # informations end here

                # allocate space for the motif probability matrix
                motif_probs = pd.DataFrame(index=alphabet, columns=range(motif_width),
                                            data=np.double(0))

                motif_width_lst.append(motif_width)
                site_counts_lst.append(site_counts)
                motif_probs_lst.append(motif_probs)

                datastart = True  # at next step begin data

                # initialize nucleotide data
                a = []
                c = []
                g = []
                t = []
                continue

            if datastart and pos_read < motif_width:
                freqs = line.split()
                a.append(np.double(freqs[0]))
                c.append(np.double(freqs[1]))
                g.append(np.double(freqs[2]))
                t.append(np.double(freqs[3]))  
                pos_read += 1

            # we read all the motif
            if pos_read == motif_width:
                a_lst.append(a)
                c_lst.append(c)
                g_lst.append(g)
                t_lst.append(t)

                # update stats about found motifs
                motifs_found += 1

                # clear the statistics
                pos_read = 0
                motif_width = None
                datastart = False
                alphalen = -1
                datastart = False 

    except:  # something went wrong
        msg = ' '.join(["Unable to read file", motif_file])
        raise FileReadingException(msg)
        die(1)

    else:

        # read the background
        if bg_file:
            bgs = readBGfile(bg_file)
        else:
            bgs = get_uniformBG(alphabet)

        bgs = pseudo_bg(bgs, no_reverse)

        motif_lst = [] # list of motifs found

        for i in range(motifs_found):
            mp = motif_probs_lst[i]

            mp.loc['A'] = a_lst[i]
            mp.loc['C'] = c_lst[i]
            mp.loc['G'] = g_lst[i]
            mp.loc['T'] = t_lst[i]

            mw = motif_width_lst[i]
            sc = site_counts_lst[i]

            mp = norm_motif(mp, mw, alphabet)
            mp = apply_pseudocount_meme(mp, pseudocount, sc, mw, bgs, alphabet)

            motif = Motif(mp, mw, alphabet, motifID_lst[i], motifName_lst[i])
            motif.setBg(bgs)

            motif_lst.append(motif)

        return motif_lst

    finally:
        mf.close() # close the file anyway


def scale_pwm(motif_matrix, alphabet, motif_width):
    """
        Scale the motif matrix values
        ----
        Parameters:
            motif_matrix (str) : path to the motif file
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

    min_val = min(motif_matrix.min())
    max_val = max(motif_matrix.max())
    motif_matrix_sc = pd.DataFrame(index=list(motif_matrix.index), columns=list(motif_matrix.columns),
                                    data=0)

    lower = min_val
    upper = max_val
        
    if lower == upper: # all values are equal
        lower = np.double(upper-1)
        
    lower = np.floor(lower)
    offset = np.round(np.floor(lower))
    scale_factor = np.floor(RANGE/(upper-lower))

    # values will be in [0, 1000]
    for nuc in alphabet:
        for j in range(motif_width):
            scaled_score = np.round((motif_matrix.loc[nuc, j]-(offset))*scale_factor)
            motif_matrix_sc.loc[nuc, j] = scaled_score
    
    # make sure the values are integers
    motif_matrix_sc[:] = motif_matrix_sc[:].astype(int)

    # now they are scaled
    min_val = min(motif_matrix_sc.min())
    max_val = max(motif_matrix_sc.max())
     
    return motif_matrix_sc, min_val, max_val, int(scale_factor), offset


def get_motif_pwm(motif_file, bgs, pseudo, no_reverse):
    """
        Build the Motif object starting from the given PWM
        ----
        Parameters:
            motif_file (str) : path to the motif file
            bgs (str) : path to the background file
            pseudo (float) : value to add to the motif counts (to avoid division 
            by 0)
            no_reverse (bool) : flag parameter to consider or not the reverse 
                                complement building the Motif object
        ----
        Returns:
            motif (Motif) : Motif object
    """

    if not motif_file:
        raise MissingFileException("The motif file is missing")
        die(1)

    if (not isMEME_ff(motif_file)) and (not isJaspar_ff(motif_file)):
        raise NotValidFFException("The motif file must be in MEME or JASPAR format")
        die(1)

    if isJaspar_ff(motif_file):
        motif = build_motif_JASPAR(motif_file, bgs, pseudo, no_reverse)
        
    elif isMEME_ff(motif_file):
        motif = build_motif_MEME(motif_file, bgs, pseudo, no_reverse)

    else:
        msg = ' '.join(["Do not know what to do with file", motif_file])
        raise NotValidFFException(msg)
        die(1)

    if not isinstance(motif, list):
        motif = [motif]

    return motif


def pseudo_bg(bgs, no_reverse):
    """ 
        Apply the pseudo count to the background frequencies
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
        bgs_avg = average_bg_with_rc(bgs)
    else:
        bgs_avg = bgs
        
    bgs_proc = norm_bg(bgs_avg)

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

    bgs_avg = {}
    
    for nuc in bgs.keys():
        rc = REV_COMPL[nuc]

        if REV_COMPL[rc] == nuc and ord(nuc) < ord(rc):
            avg_freq = np.double((bgs[nuc]+bgs[rc])/np.double(2))
            bgs_avg.update({nuc: avg_freq})
            bgs_avg.update({rc:avg_freq})
            
    return bgs_avg


def norm_bg(bgs):
    """
        Normalize the background frequencies
        ----
        Parameters:
            Parameters:
            bgs (dict) : dictionary of the background probabilities
        ----
        Returns:
            bgs_norm (dict) : normalized background probabilities
    """
    
    PSEUDO = np.double(0.0000005) # pseudocount

    alphabet = sorted(list(bgs.keys()))
    tot = np.double(len(alphabet)*PSEUDO)
    bgs_norm = {}
    
    for nuc in bgs.keys():
        tot += np.double(bgs[nuc])

    assert tot > 0
        
    for nuc in bgs.keys():
        prob = np.double((bgs[nuc]+PSEUDO)/tot)
        bgs_norm.update({nuc : prob})
        
    tot = np.double(0)
    for nuc in bgs.keys():
        tot += bgs[nuc]
    
    assert tot != 0
            
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
    
    tolerance = 0.00001 # allowed tolerance in the difference between the position
                        # probability and 1
                      
    for j in range(motif_width):
        tot = np.double(0)

        for nuc in alphabet:
            tot += motif_probs.loc[nuc, j]

        assert tot != 0
            
        if not almost_equal(1, tot, tolerance):
            for nuc in alphabet:
                motif_probs.loc[nuc, j] = np.double(motif_probs.loc[nuc, j]/tot)
                
    return motif_probs

