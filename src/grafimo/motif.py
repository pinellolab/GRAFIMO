"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

The file contains the definition of the Motif class, that represent
the motif that will be searched on the genome.

"""


from grafimo.GRAFIMOException import NoDataFrameException, WrongMotifWidthException, WrongMotifIDException, \
        WrongMotifNameException, NotValidMotifMatrixException, NotValidBGException, \
        NotValidAlphabetException, NotValidFFException, FileReadingException, \
        ValueException
from grafimo.utils import die, DNA_ALPHABET, REV_COMPL, PSEUDObg, isListEqual, isJaspar_ff, isMEME_ff, \
        almost_equal, RANGE, printProgressBar, sigint_handler
from motif_processing import readBGfile, get_uniformBG, apply_pseudocount_jaspar, apply_pseudocount_meme, \
        compute_log_odds, comp_pval_mat
import multiprocessing as mp
import signal
import pandas as pd
import numpy as np
import time
import os


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

    _count_matrix = None  # pd.DataFrame
    _score_matrix = None  # np.ndarray
    _pval_matrix = None  # np.array
    _min_val = -np.inf  # np.double
    _max_val = np.inf  # np.double
    _scale = -1  # int
    _offset = 0  # np.double
    _bg = None  # background distribution (dict)
    _width = -1  # int
    _motif_id = None  # jaspar ID
    _motif_name = None  # common TF name
    _alphabet = None
    _isScaled = False

    def __init__(self,
                 count_matrix,
                 width,
                 alphabet,
                 motif_id,
                 motif_name):

        if count_matrix.empty:
            errmsg = "\n\nERROR: attempt to initialize the motif object with an empty count matrix"
            raise NotValidMotifMatrixException(errmsg)

        if not isinstance(count_matrix, pd.DataFrame):
            raise NoDataFrameException("\n\nERROR: the given value is not a pandas.DatFrame instance")

        if not isinstance(width, int) or width < 0:
            errmsg = "\n\nERROR: attempt to initialize motif without a valid width"
            raise WrongMotifWidthException(errmsg)

        if not isinstance(motif_id, str) or not motif_id:
            raise WrongMotifIDException("\n\nERROR: cannot initialize the motif with the given ID")

        if not isinstance(motif_name, str) or not motif_name:
            raise WrongMotifNameException("\n\nERROR: cannot initialize the motif with the given name")

        if not isinstance(alphabet, list) or not isListEqual(alphabet, DNA_ALPHABET):
            errmsg = "\n\nERROR: cannot initialize a motif object with a wrong alphabet"
            raise NotValidAlphabetException(errmsg)

        self._count_matrix = count_matrix
        self._width = width
        self._motif_id = motif_id
        self._motif_name = motif_name
        self._alphabet = alphabet
    # end of __init__()

    ### setter methods ###

    def setMotif_matrix(self,
                        motif_matrix):

        if motif_matrix.empty:
            errmsg = "\n\nERROR: attempt to use an empty motif matrix"
            raise NotValidMotifMatrixException(errmsg)

        if not isinstance(motif_matrix, pd.DataFrame):
            raise NoDataFrameException("\n\nERROR: the given value is not a pandas.DataFrame instance")

        self._count_matrix = motif_matrix

    def setMotif_scoreMatrix(self,
                             score_matrix):

        if (not isinstance(score_matrix, np.ndarray) and
                not isinstance(score_matrix, pd.DataFrame)):
            errmsg = "\n\nERROR: the given data-structure is not an instance of numpy.ndarray or pandas.DataFrame"
            raise ValueError(errmsg)

        if isinstance(score_matrix, pd.DataFrame):
            if score_matrix.empty:
                errmsg = "\n\nERROR: attempt to use an empty score matrix"
                raise NotValidMotifMatrixException(errmsg)

        if isinstance(score_matrix, np.ndarray):
            if score_matrix.size == 0:
                errmsg = "\n\nERROR: attempt to use an empty score matrix"
                raise NotValidMotifMatrixException(errmsg)

        self._score_matrix = score_matrix

    def setMotif_pval_matrix(self,
                             pval_mat):

        if len(pval_mat) == 0 or sum(pval_mat[:]) <= 0:  # empty or not valid p-value matrix
            raise NotValidMotifMatrixException("\n\nERROR: invalid p-value matrix")

        self._pval_matrix = pval_mat

    def setMin_val(self,
                   min_val):

        if min_val <= -np.inf:
            errmsg = ' '.join(["\n\nERROR: impossible to assign", min_val, "to Motif.min_val"])
            raise ValueException(errmsg)

        self._min_val = min_val

    def setMax_val(self,
                   max_val):

        if max_val >= np.inf:
            errmsg = ' '.join(["\n\nERROR: impossible to assign", max_val, "to Motif.max_val"])
            raise ValueException(errmsg)

        self._max_val = max_val

    def setScale(self,
                 scale):

        if not isinstance(scale, int):
            raise ValueException("\n\nERROR: the scale factor must be an int")

        assert scale > 0

        self._scale = scale

    def setOffset(self,
                  offset):

        self._offset = offset

    def setBg(self,
              bgs):

        if not isinstance(bgs, dict):
            raise NotValidBGException("\n\nERROR: the background values are not in a dictionary")

        self._bg = bgs

    def setWidth(self,
                 width):

        if not isinstance(width, int) or width <= 0:
            errmsg = "\n\nERROR: attempt to initialize motif without a valid width"
            raise WrongMotifWidthException(errmsg)

        self._width = width

    def setMotifID(self,
                   motif_id):

        if not isinstance(motif_id, str):
            errmsg = ' '.join(["\n\nERROR: cannot initialize the motif with the given ID:", motif_id])
            raise WrongMotifIDException(errmsg)

        if not motif_id:
            raise WrongMotifIDException("\n\nERROR: cannot use a motif with an empty ID")

        self._motif_id = motif_id

    def setMotifName(self,
                     motif_name):

        if not isinstance(motif_name, str):
            errmsg = ' '.join(["Cannot initialize the motif with the given name:", motif_name])
            raise WrongMotifNameException(errmsg)

        if not motif_name:
            raise WrongMotifNameException("\n\nERROR: cannot use a motif with an empty name")

        self._motif_name = motif_name

    def setAlphabet(self,
                    alphabet):

        if not isinstance(alphabet, list):
            raise NotValidAlphabetException("\n\nERROR: the given alphabet is not in a list")

        if not isListEqual(alphabet, DNA_ALPHABET):
            raise NotValidAlphabetException("\n\nERROR: the given alphabet is not a valid DNA alphabet")

        self.alphabet = alphabet

    def setIsScaled(self,
                    isScaled):

        if not isinstance(isScaled, bool):
            raise Exception("\n\nERROR: the isScaled value must be a boolean")

        self._isScaled = isScaled

    ### getter methods ###

    def getMotif_matrix(self):

        return self._count_matrix

    def getMotif_scoreMatrix(self):

        return self._score_matrix

    def getMotif_pval_mat(self):

        return self._pval_matrix

    def getMin_val(self):

        return self._min_val

    def getMax_val(self):

        return self._max_val

    def getScale(self):

        return self._scale

    def getOffset(self):

        return self._offset

    def getBg(self):

        return self._bg

    def getWidth(self):

        return self._width

    def getMotifID(self):

        return self._motif_id

    def getMotifName(self):

        return self._motif_name

    def getAlphabet(self):

        return self._alphabet

    def getIsScaled(self):

        return self._isScaled

    def compute_minValue(self):

        motif_matrix = self.getMotif_matrix()
        min_value = motif_matrix.min().sum()
        self._min_val = min_value

    def print(self,
              matrix):
        allowed_matrices = ["raw_counts", "score_matrix", "pval_matrix"]

        if matrix not in allowed_matrices:
            raise ValueError("ERROR: unknown Motif matrix to print")

        if str(matrix) == "raw_counts":
            print(self._count_matrix)
        elif str(matrix) == "score_matrix":
            print(self._score_matrix)
        elif str(matrix) == "pval_matrix":
            print(self._pval_matrix)
        else:
            # we should not reach this point
            raise ValueError("ERROR: unknown Motif matrix to print")
    # end of print()
# end of Motif


#################################################
# end of Motif definition
#################################################


def build_motif_JASPAR(motif_file,
                       bg_file,
                       pseudocount,
                       no_reverse,
                       verbose):
    """
        Build a Motif object starting from raw counts
        data stored in a JASPAR motif file.

        The raw counts are processed and the resulting values
        are used to define the scoring matrix for the motif
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

    if not motif_file:
        raise FileNotFoundError("\n\nERROR: the motif file is missing")
        die(1)

    # check if the input file is in JASPAR format
    if not isJaspar_ff(motif_file):
        raise NotValidFFException("ERROR: the given motif file is not in JASPAR or MEME format")
        die(1)

    assert pseudocount > 0

    # read the motif file
    motif = read_JASPAR_motif(motif_file, bg_file, pseudocount, no_reverse, verbose)

    if verbose:
        start_mp = time.time()

    #  get log-odds values for motif
    motif = process_motif_for_logodds(motif)

    if verbose:
        end_mp = time.time()
        msg = ''.join(["Processed motif ", motif.getMotifID(), " in ", str(end_mp - start_mp), "s"])
        print(msg)
    # end if

    return motif
# end of build_motif_JASPAR()


def read_JASPAR_motif(motif_file,
                      bg_file,
                      pseudocount,
                      no_reverse,
                      verbose):
    """
        Read data contained in a JASPAR motif file and build a Motif
        object from them
        ----
        Params:
            motif_file (str) : path to the motif file (in JASPAR format)
            bg_file (str) : path to the background file
            no_reverse (bool) : flag parameter to consider or not the reverse
                                complement building the Motif object
        ----
        Returns:
            motif (Motif) : Motif object summarizing data contained in
                                motif_file
    """

    # lists where store nucleotides and raw counts
    nucs = []
    counts = []

    if verbose:
        start_rm = time.time()

    try:
        # open the motif file
        with open(motif_file) as in_mtf:

            header = str(in_mtf.readline()[1:])  # read the header
            motifID, motifName = header.split('\t')[0:2]  # get the jaspar ID and the common TF name
            motifName = motifName[:-1]  # remove '\n'

            for line in in_mtf:
                line = line.strip()
                nuc = line.strip()[:1]  # read nucleotide
                count = list(map(float, line.strip()[1:].split()[1:][:-1]))  # read raw counts

                nucs.append(nuc)
                counts.append(count)
            # end for
        # end open

    except:  # something went wrong
        errmsg = ' '.join(["\n\nERROR: unable to read file", motif_file])
        raise FileReadingException(errmsg)

    else:

        motif_counts = pd.DataFrame(data=counts, index=nucs)  # raw counts
        motif_width = int(len(counts[0]))  # the check of equal length for all raw counts is made building the DataFrame
        alphabet = sorted(nucs)  # alphabet as list

        # read the background file
        if bg_file == 'UNIF':
            bgs = get_uniformBG(alphabet)
        elif os.path.exists(bg_file):
            bgs = readBGfile(bg_file)
        else:
            raise NotValidBGException("\n\nERROR: unable to find the given background file")
        # end if

        bgs = pseudo_bg(bgs, no_reverse)

        motif_probs = (motif_counts / motif_counts.sum(0))  # get probabilities
        motif_probs = norm_motif(motif_probs, motif_width, alphabet)
        motif_probs = apply_pseudocount_jaspar(motif_counts, motif_probs, pseudocount, bgs,
                                               motif_width, alphabet)

        motif = Motif(motif_probs, motif_width, alphabet, motifID, motifName)
        motif.setBg(bgs)

        if verbose:
            end_rm = time.time()
            msg = ''.join(["Read motif ", motifID, " in ", str(end_rm - start_rm), "s"])
            print(msg)
        # end if

        return motif

    finally:
        in_mtf.close()  # close the motif file anyway
# end of read_JASPAR_motif()


def build_motif_MEME(motif_file,
                     bg_file,
                     pseudocount,
                     no_reverse,
                     cores,
                     verbose):
    """
        Build a the Motif object starting from the data
        stored in a given MEME file.

        The probabilities are processed and the resulting values
        are used to build the scoring matrix for the motif.
        ----
        Parameters:
            motif_file (str) : path to the motif file
            bg_file (str) : path to the background file
            pseudocount (float) : value to add to the motif counts
            no_reverse (bool) : if set to True, only data related to
                                forward strand will be used
            cores (int) : number of cores to use, during motif processing
        ----
        Returns:
            motif (Motif) : Motif object built from data contained in
                            motif_file
    """

    if not motif_file:
        raise FileNotFoundError("\n\nERROR: the motif file is missing")

    # check if the input is in MEME format
    if not isMEME_ff(motif_file):
        # if in other format we should not be here
        raise NotValidFFException("\n\nERROR: the given motif file is not in MEME format")

    if verbose:
        start_rm_all = time.time()

    # read the motif file
    motif_lst = read_MEME_motif(motif_file, bg_file, pseudocount, no_reverse, verbose)
    motif_num = len(motif_lst)

    if verbose:
        end_rm_all = time.time()
        msg = ''.join(["\nRead all motif contained in ", motif_file, " in ", str(end_rm_all - start_rm_all), "s"])
        print(msg)
    # end if

    print("\nRead", motif_num, "motifs in", motif_file)
    print("\nProcessing motifs\n")

    # list of the fully processed motifs
    complete_motifs = []

    if verbose:
        start_mp_all = time.time()

    # process each found motif
    if motif_num >= cores:  # worth to use multiprocessing

        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool = mp.Pool(processes=cores)  # use #cores processes
        signal.signal(signal.SIGINT, original_sigint_handler)  # overwrite the default SIGINT handler to exit gracefully
        # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python

        try:
            res = (pool.map_async(process_motif_for_logodds, motif_lst))

            it = 0
            while (True):
                if res.ready():
                    # when finished call for the last time printProgressBar()
                    printProgressBar(tot, tot, prefix='Progress:',
                                        suffix='Complete', length=50)
                    break
                # end if
                if it == 0:
                    tot = res._number_left

                remaining = res._number_left
                printProgressBar((tot - remaining), tot, prefix='Progress:',
                                suffix='Complete', length=50)
                time.sleep(2)
                it += 1
            # end while

            complete_motifs += res.get(60 * 60 * 60)  # does not ignore signals

        except KeyboardInterrupt:
            pool.terminate()
            sigint_handler()

        else:
            pool.close()

            if verbose:
                end_mp_all = time.time()
                msg = ''.join(["Processed all motifs contained in ", motif_file, " in ",
                               str(end_mp_all - start_mp_all), "s"])
                print(msg)
            # end if

            return complete_motifs
        # end try

    else:  # the sequential execution is fine

        # process each found motif
        for m in motif_lst:
            complete_motifs.append(process_motif_for_logodds(m))

        if verbose:
            end_mp_all = time.time()
            msg = ''.join(["Processed all motifs contained in ", motif_file, " in ",
                           str(end_mp_all - start_mp_all), "s"])
            print(msg)
        # end if

        return complete_motifs

# end build_motif_MEME()


def read_MEME_motif(motif_file,
                    bg_file,
                    pseudocount,
                    no_reverse,
                    verbose):
    """
        Read the motif file in MEME format and build a motif
        object from it.
        Note that a MEME file can contain a variable number of
        motifs
        ----
        Params:
            motif_file (str) : path to the motif file
            bg_file (str) : path to the background file
            pseudocount (np.double) : pseudocount to add to motif frequencies
            no_reverse (bool) : if set to True, only data related to
                                forward strand will be used
        ----
        Returns:
            motif (Motif) : returns a Motif object
    """

    try:
        with open(motif_file, 'r') as in_mtf:  # open the motif file

            infostart = False  # flag to keep track were the infos about the motif begin
            datastart = False  # flag to keep track were the motif data begin
            motifs_found = 0  # number of motifs found in the MEME file

            motifID_lst = []  # list of the found motif IDs
            motifName_lst = []  # list of the found motif names
            motif_width_lst = []  # list of the found motif widths
            site_counts_lst = []  # list of the found motif site counts
            alphalen_lst = []  # list of the found motif alphabet lengths
            motif_probs_lst = []  # list of the found motif probability matrices
            a_lst = []  # list of the found As probabilities for each motif
            c_lst = []  # list of the found Cs probabilities for each motif
            g_lst = []  # list of the found Gs probabilities for each motif
            t_lst = []  # list of the found Ts probabilities for each motif

            motif_width = None
            pos_read = 0

            for line in in_mtf:
                if line[0:8] == 'ALPHABET':
                    alphabet = sorted(list(set(line[10:-1])))
                    assert isListEqual(alphabet, DNA_ALPHABET)

                if line[0:5] == 'MOTIF':

                    if verbose:
                        start_rm = time.time()

                    motifID, motifName = line.split()[1:3]

                    motifID_lst.append(motifID)
                    motifName_lst.append(motifName)

                    # the informations about motif start here
                    infostart = True
                    continue
                # end if

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
                # end if

                if datastart and pos_read < motif_width:
                    freqs = line.split()
                    a.append(np.double(freqs[0]))
                    c.append(np.double(freqs[1]))
                    g.append(np.double(freqs[2]))
                    t.append(np.double(freqs[3]))
                    pos_read += 1
                # end if

                # we read all current motif data
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

                    if verbose:
                        end_rm = time.time()
                        msg = ''.join(["Read motif ", motifID, " in ", str(end_rm - start_rm), "s"])
                        print(msg)
                    # end if
                # end if

    except:  # something went wrong
        errmsg = ' '.join(["Unable to read file", motif_file])
        raise FileReadingException(errmsg)

    else:

        # read the background
        if bg_file == 'UNIF':
            bgs = get_uniformBG(alphabet)
        elif os.path.exists(bg_file):
            bgs = readBGfile(bg_file)
        else:
            raise NotValidBGException("\n\nERROR: unable to find the given background file")
        # end if

        bgs = pseudo_bg(bgs, no_reverse)

        motif_lst = []  # list of found motifs

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
        # end for

        return motif_lst

    finally:
        in_mtf.close()  # close the file anyway
    # end try
# end read_MEME_motif()


def process_motif_for_logodds(motif):
    """
        Process the given motif and return the
        logodds values from motif file data.

        It is also computed the P-value matrix, using
        a DP-algorithm (Staden, R., 1994)
        ----
        Parameters:
            motif (Motif): motif to process
        ----
        Return:
            motif (Motif): the update input motif
    """

    # get the log-odds
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

    motif.setMotif_scoreMatrix(scaled_scores.values)

    return motif
# end of process_motif_for_logodds()


def scale_pwm(motif_matrix,
              alphabet,
              motif_width):
    """
        Scale the motif matrix values
        ----
        Parameters:
            motif_matrix (str) : count matrix
            alphabet (str) : motif alphabet
            motif_width (int) : motif width
        ----
        Returns:
            motif_matrix_sc (np.ndarray) : scaled motif matrix
            min_val (int) : lowest value in the scaled motif matrix
            max_val (int) : higest value in the scaled motif matrix
            scale_factor (int)
            offset (int)
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

    if lower == upper:  # all values are equal
        lower = np.double(upper - 1)

    lower = np.floor(lower)
    offset = np.round(np.floor(lower))
    scale_factor = np.floor(RANGE / (upper - lower))

    # values will be in [0, 1000]
    for nuc in alphabet:
        for j in range(motif_width):
            scaled_score = np.round((motif_matrix.loc[nuc, j] - (offset)) * scale_factor)
            motif_matrix_sc.loc[nuc, j] = scaled_score
        # end for
    # end for

    # make sure the values are integers
    motif_matrix_sc[:] = motif_matrix_sc[:].astype(int)

    # now they are scaled
    min_val = min(motif_matrix_sc.min())
    max_val = max(motif_matrix_sc.max())

    return motif_matrix_sc, min_val, max_val, int(scale_factor), offset
# end of scale_pwm()


def get_motif_pwm(motif_file,
                  args_obj,
                  cores):
    """
        Build a Motif object starting from a given PWM
        ----
        Parameters:
            motif_file (str) : motif file to process
            args_obj (Findmotif): data-structure containing the
                                    parameters to scan a given
                                    VG or a set of VGs
            cores (int) : number of cores to use during motif
                            processing
        ----
        Returns:
            motif (list) : list of processed motifs as Motif objects
    """

    # get arguments required to process the motif
    bgs = args_obj.get_bgfile()
    pseudo = args_obj.get_pseudo()
    no_reverse = args_obj.get_no_reverse()
    verbose = args_obj.get_verbose()

    if not motif_file:
        raise FileNotFoundError("\n\nERROR: the motif file is missing")

    if (not isMEME_ff(motif_file)) and (not isJaspar_ff(motif_file)):
        raise NotValidFFException("\n\nERROR: the motif file must be in MEME or JASPAR format")

    if isJaspar_ff(motif_file):
        motif = build_motif_JASPAR(motif_file, bgs, pseudo, no_reverse, verbose)

    elif isMEME_ff(motif_file):
        motif = build_motif_MEME(motif_file, bgs, pseudo, no_reverse, cores, verbose)

    else:
        errmsg = ' '.join(["\n\nERROR: do not know what to do with file", motif_file])
        raise NotValidFFException(errmsg)
    # end if

    if not isinstance(motif, list):
        motif = [motif]

    return motif
# end of get_motif_pwm()


def pseudo_bg(bgs,
              no_reverse):
    """
        Add the pseudocount to the background frequencies
        ----
        Parameters:
            bgs (dict) : dictionary of the background frequencies
            no_reverse (bool) : if set to True, the background
                                frequencies will be averaged with the
                                reverse complement frequencies
        ----
        Returns:
            bgs_proc (dict) : normalized (and averaged) background frequencies
    """

    if not isinstance(bgs, dict):
        raise NotValidBGException("\n\nERROR: unable to add the pseudocount to the background")

    if not isinstance(no_reverse, bool):
        raise ValueException(' '.join(["Boolean value required, got", str(type(no_reverse))]))

    if not no_reverse:
        bgs_avg = average_bg_with_rc(bgs)
    else:
        bgs_avg = bgs

    bgs_proc = norm_bg(bgs_avg)

    return bgs_proc
# end pseudo_bg()


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
            avg_freq = np.double((bgs[nuc] + bgs[rc]) / np.double(2))
            bgs_avg.update({nuc: avg_freq})
            bgs_avg.update({rc: avg_freq})
        # end if
    # end for

    return bgs_avg
# end average_bg_with_rc()


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

    # PSEUDO = np.double(0.0000005)  # pseudocount

    alphabet = sorted(list(bgs.keys()))
    tot = np.double(len(alphabet) * PSEUDObg)
    bgs_norm = {}

    for nuc in bgs.keys():
        tot += np.double(bgs[nuc])

    assert tot > 0

    for nuc in bgs.keys():
        prob = np.double((bgs[nuc] + PSEUDObg) / tot)
        bgs_norm.update({nuc: prob})

    tot = np.double(0)
    for nuc in bgs.keys():
        tot += bgs[nuc]

    assert tot != 0

    return bgs_norm
# end norm_bg()


def norm_motif(motif_probs,
               motif_width,
               alphabet):
    """
        Normalize motif probabilities
        ----
        Parameters:
            motif_probs (pd.DataFarme) : probability matrix
            motif_width (int) : motif width
            alphabet (list) : motif alphabet
        ----
        Returns:
            motif_probs (pd.DataFrame) : normalized probability matrix
    """

    tolerance = 0.00001  # allowed tolerance in the difference between the position probability and 1

    for j in range(motif_width):
        tot = np.double(0)

        for nuc in alphabet:
            tot += motif_probs.loc[nuc, j]

        assert tot != 0

        if not almost_equal(1, tot, tolerance):
            for nuc in alphabet:
                motif_probs.loc[nuc, j] = np.double(motif_probs.loc[nuc, j] / tot)
        # end if
    # end for

    return motif_probs
# end norm_motif()

