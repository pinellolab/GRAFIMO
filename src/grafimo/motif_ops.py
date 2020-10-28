"""Functions to process the input motif Position Weight Matrix.

The PWM is processed in order to get log-odds values for each entry and
have a log-likelihood score to assign to the motif occurrences 
candidates.

The resulting scoring matrix is then assigned to a Motif object, which 
contains all the information useful to describe and summarize the given 
motif PWM

"""


from grafimo.GRAFIMOException import NoDataFrameException, ValueException, \
    NotValidMotifMatrixException, NotValidBGException, FileReadingException, \
    NotValidAlphabetException, NotValidFFException
from grafimo.utils import die, DNA_ALPHABET, REV_COMPL, PSEUDObg, isListEqual,\
     isJaspar_ff, isMEME_ff, almost_equal, RANGE, printProgressBar, \
     sigint_handler
from motif_processing import readBGfile, get_uniformBG, \
    apply_pseudocount_jaspar, apply_pseudocount_meme, compute_log_odds, \
    comp_pval_mat
from grafimo.workflow import Findmotif
from typing import List, Optional, Dict, Tuple
from grafimo.motif import Motif
import multiprocessing as mp
import pandas as pd
import numpy as np
import signal
import time
import os


def build_motif_JASPAR(motif_file: str,
                       bg_file: str,
                       pseudocount: float,
                       no_reverse: bool,
                       verbose: bool
) -> Motif:
    """Build the Motif object from a JASPAR motif Position Weight
    Matrix.

    Is computed the scoring matrix from the values given with the PWM
    and the P-value matrix to assign a statistical significance to
    each motif occurrence candidate, based on the resulting log-odds
    score.

    Parameters
    ----------
    motif_file : str
        path to the motif PWM
    bg_file
        path to the background file in Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html).
    pseudocount : float
        value to add to motif PWM counts
    no_reverse : bool
        if False only the forward strand will be considered, otherwise
        both forward and reverse are considered
    verbose : bool
        print additional information

    Returns
    -------
    Motif
        processed motif object
    """

    errmsg: str
    if not motif_file:
        errmsg = "\n\nERROR: the motif file is missing"
        raise FileNotFoundError(errmsg)

    if not isJaspar_ff(motif_file):
        errmsg = "ERROR: the given motif file is not in JASPAR or MEME format"
        raise NotValidFFException(errmsg)

    assert pseudocount > 0

    # read the motif file
    motif: Motif = read_JASPAR_motif(motif_file, bg_file, pseudocount, 
                                    no_reverse, verbose)

    if verbose:
        start_mp: float = time.time()

    #  get log-odds values for motif
    motif = process_motif_for_logodds(motif)

    if verbose:
        end_mp: float = time.time()
        print(
            "Processed motif %s in %.2fs" % (motif.getMotifID(), 
                                              (end_mp - start_mp))
             )
    # end if

    return motif

# end of build_motif_JASPAR()


def read_JASPAR_motif(motif_file: str,
                      bg_file: str,
                      pseudocount: float,
                      no_reverse: bool,
                      verbose: bool
) -> Motif:
    """Read a motif PWM in JASPAR format.

    The data read are then used to build the scoring matrix for the 
    motif, the P-value matrix, etc.

    Parameters:
    motif_file : str
        path to the motif PWM in JASPAR format
    bg_file
        path to the background file in Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html).
    pseudocount : float
        value to add to motif PWM counts
    no_reverse : bool
        if False only the forward strand will be considered, otherwise
        both forward and reverse are considered
    verbose : bool
        print additional information

    Returns
    -------
    Motif
        Motif object storing the data contained in motif_file
    """

    nucs: List[str]
    counts: List[float]

    # lists where store nucleotides and raw counts
    nucs = list()
    counts = list()

    if verbose:
        start_rm: float = time.time()

    try:
        # open the motif file
        with open(motif_file) as in_mtf:

            header: str
            motifID: str
            motifName: str
            
            # read the header
            header = str(in_mtf.readline()[1:])  
            # get the jaspar ID and the common TF name
            motifID, motifName = header.split('\t')[0:2]  
            motifName = motifName[:-1]  # remove '\n'

            for line in in_mtf:
                line = line.strip()
                nuc = line.strip()[:1] 
                count = list(map(float, line.strip()[1:].split()[1:][:-1]))  

                nucs.append(nuc)
                counts.append(count)
            # end for
        # end open

    except:
        errmsg: str = ' '.join(["\n\nERROR: unable to read file", motif_file])
        raise FileReadingException(errmsg)

    else:

        motif_counts = pd.DataFrame(data=counts, index=nucs)  
        # the check of equal length for all raw counts is made building 
        # the DataFrame
        motif_width: int = int(len(counts[0]))  
        alphabet: list = sorted(nucs)  # alphabet as list

        bgs: Dict
        if bg_file == 'UNIF':
            bgs = get_uniformBG(alphabet)
        elif os.path.exists(bg_file):
            bgs = readBGfile(bg_file)
        else:
            errmsg = "\n\nERROR: unable to find the given background file"
            raise NotValidBGException(errmsg)
        # end if

        bgs = pseudo_bg(bgs, no_reverse)

        motif_probs: pd.DataFrame
        motif_probs = (motif_counts / motif_counts.sum(0)) 
        motif_probs = norm_motif(motif_probs, motif_width, alphabet)
        motif_probs = apply_pseudocount_jaspar(motif_counts, motif_probs, 
                                               pseudocount, bgs, motif_width, 
                                               alphabet)

        motif: Motif = Motif(motif_probs, motif_width, alphabet, motifID,
                             motifName)
        motif.setBg(bgs)

        if verbose:
            end_rm: float = time.time()
            msg: str = ''.join(["Read motif ", motifID, " in ", 
                                str(end_rm - start_rm), "s"])
            print(msg)
        # end if

        return motif

    finally:
        in_mtf.close()  # close the motif file anyway

# end of read_JASPAR_motif()


def build_motif_MEME(motif_file: str,
                     bg_file: str,
                     pseudocount: float,
                     no_reverse: bool,
                     cores: int,
                     verbose: bool
) -> List[Motif]:
    """Read a motif PWM in MEME format.

    The data read are then used to build the scoring matrix for the 
    motif, the P-value matrix, etc.

    Parameters:
    motif_file : str
        path to the motif PWM in MEME format
    bg_file
        path to the background file in Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html).
    pseudocount : float
        value to add to motif PWM counts
    no_reverse : bool
        if False only the forward strand will be considered, otherwise
        both forward and reverse are considered
    cores : int
        Number of cores to use while building the Motif object
    verbose : bool
        print additional information

    Returns
    -------
    Motif
        Motif object storing the data contained in motif_file
    """

    errmsg: str
    if not motif_file:
        errmsg = "\n\nERROR: the motif file is missing"
        raise FileNotFoundError(errmsg)

    if not isMEME_ff(motif_file):
        errmsg = "\n\nERROR: the given motif file is not in MEME format"
        raise NotValidFFException(errmsg)

    if verbose:
        start_rm_all: float = time.time()

    motif_lst: List[Motif]
    motif_lst = read_MEME_motif(motif_file, bg_file, pseudocount, no_reverse, 
                                verbose)
    motif_num: int = len(motif_lst)

    if verbose:
        end_rm_all: float = time.time()
        msg: str = ''.join(["\nRead all motif contained in ", motif_file, 
                            " in ", str(end_rm_all - start_rm_all), "s"])
        print(msg)
    # end if

    print("\nRead", motif_num, "motifs in", motif_file)
    print("\nProcessing motifs\n")

    # list of the fully processed motifs
    complete_motifs = list()

    if verbose:
        start_mp_all: str = time.time()

    # process each found motif
    if motif_num >= cores:  # worth to use multiprocessing

        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool: mp.Pool = mp.Pool(processes=cores)  # use #cores processes
        signal.signal(signal.SIGINT, original_sigint_handler)  
        # overwrite the default SIGINT handler to exit gracefully
        # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python

        try:
            res = (pool.map_async(process_motif_for_logodds, motif_lst))

            it: int = 0
            while (True):
                if res.ready():
                    # when finished call for the last time 
                    # printProgressBar()
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

            # does not ignore signals
            complete_motifs += res.get(60 * 60 * 60)  

        except KeyboardInterrupt:
            pool.terminate()
            sigint_handler()

        else:
            pool.close()

            if verbose:
                end_mp_all: float = time.time()
                print(
                    "Processed motif(s) in %s in %.2fs" % 
                    (motif_file, (end_mp_all - start_mp_all))
                    )
            # end if

            return complete_motifs
        # end try

    else: 

        # process each found motif
        for m in motif_lst:
            complete_motifs.append(process_motif_for_logodds(m))

        if verbose:
            end_mp_all: float = time.time()
            print(
                "Processed motif(s) in %s in %.2fs" % 
                (motif_file, (end_mp_all - start_mp_all))
                )
        # end if

        return complete_motifs

# end build_motif_MEME()


def read_MEME_motif(motif_file: str,
                    bg_file: str,
                    pseudocount: float,
                    no_reverse: bool,
                    verbose: bool
) -> List[Motif]:
    """Read a motif PWM in MEME format.

    The data read are then used to build the scoring matrix for the 
    motif, the P-value matrix, etc.

    Since a MEME file can contain one or more motifs, for each PWM
    contained is built the corresponding motif object.
    The resulting set of motifs are then stored in a list.
    
    Parameters
    ----------
    motif_file : str
        path to the motif PWM
    bg_file : str
        path to the background probability distribution
    pseudocount : float
        pseudocount to add to the PWM values
    no_reverse : bool
        if False only the forward strand will be considered, otherwise
        both forward and reverse are considered
    verbose : bool
        print additional information

    Returns
    -------
    List[Motif]
        List of Motif objects storing the data contained in motif_file
    """

    try:
        with open(motif_file, 'r') as in_mtf:  # open the motif file
            
            # flag to keep track were the infos about the motif begin
            infostart: bool
            # flag to keep track were the motif data begin
            datastart: bool
            # number of motifs found in the MEME file
            motifs_found: int
            # list of the found motif IDs
            motifID_lst: List[str]
            # list of the found motif names
            motifName_lst: List[str]
            # list of the found motif widths
            motif_width_lst: List[int]
            # list of the found motif site counts
            site_counts_lst: List[int]
            # list of the found motif alphabet lengths
            alphalen_lst: List[int]
            # list of the found motif probability matrices
            motif_probs_lst: List[pd.DataFrame]
            # list of the found As probabilities for each motif
            a_lst: List[np.double]
            # list of the found Cs probabilities for each motif
            c_lst: List[np.double]
            # list of the found Gs probabilities for each motif
            g_lst: List[np.double]
            # list of the found Ts probabilities for each motif
            t_lst: List[np.double]

            infostart = False  
            datastart = False  
            motifs_found = 0  

            motifID_lst = list()  
            motifName_lst = list()  
            motif_width_lst = list()  
            site_counts_lst = list()  
            alphalen_lst = list()  
            motif_probs_lst = list() 
            a_lst = list()  
            c_lst = list()  
            g_lst = list() 
            t_lst = list()  
            motif_width = None
            pos_read = 0

            for line in in_mtf:
                if line[0:8] == 'ALPHABET':
                    alphabet: List = sorted(list(set(line[10:-1])))
                    assert isListEqual(alphabet, DNA_ALPHABET)

                if line[0:5] == 'MOTIF':

                    if verbose:
                        start_rm: float = time.time()

                    # read motif ID and full name
                    motif_header: str = line.split()  

                    assert len(motif_header) > 0

                    # there are two ways to define the motif name line 
                    # in MEME file
                    # (refer to http://meme-suite.org/doc/meme-format.html?man_type=web):
                    #   1 - MOTIF motif_alternate_name
                    #   2 - MOTIF motif_identifier motif_alternate_name

                    motifID: str
                    motifName: str

                    if len(motif_header) == 2:  # support case (1)
                        motifID = motif_header[1]
                        motifName = motif_header[1]

                    else:  # support case (2)
                        motifID, motifName = motif_header[1:3]
                    # end if

                    motifID_lst.append(motifID)
                    motifName_lst.append(motifName)

                    # the informations about motif start here
                    infostart = True
                    continue
                # end if

                if infostart and len(line.strip()) != 0:
                    infos: str = line[26:]
                    infosplit: List[str] = infos.split()
                    alphalen: int = int(infosplit[1])
                    alphalen_lst.append(alphalen)

                    assert alphalen == len(alphabet)

                    motif_width: int = int(infosplit[3])
                    site_counts: int = int(infosplit[5])
                    infostart = False  # informations end here

                    # allocate space for the motif probability matrix
                    motif_probs: pd.DataFrame = pd.DataFrame(index=alphabet, 
                                                    columns=range(motif_width),
                                                    data=np.double(0)
                                                            )

                    motif_width_lst.append(motif_width)
                    site_counts_lst.append(site_counts)
                    motif_probs_lst.append(motif_probs)

                    datastart = True  # at next step begin data

                    # initialize nucleotide data
                    a = list()
                    c = list()
                    g = list()
                    t = list()
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
                    pos_read: int = 0
                    motif_width = None
                    datastart = False
                    alphalen = -1
                    datastart = False

                    if verbose:
                        end_rm: float = time.time()
                        msg: str = ''.join(["Read motif ", motifID, " in ", 
                                             str(end_rm - start_rm), "s"])
                        print(msg)
                    # end if
                # end if

    except:  # something went wrong
        errmsg: str = ' '.join(["Unable to read file", motif_file])
        raise FileReadingException(errmsg)

    else:

        bgs: dict
        # read the background
        if bg_file == 'UNIF':
            bgs = get_uniformBG(alphabet)
        elif os.path.exists(bg_file):
            bgs = readBGfile(bg_file)
        else:
            errmsg = "\n\nERROR: unable to find the given background file"
            raise NotValidBGException(errmsg)
        # end if

        bgs = pseudo_bg(bgs, no_reverse)

        motif_lst: List[Motif]
        motif_lst = list()

        for i in range(motifs_found):
            mp: pd.DataFrame = motif_probs_lst[i]

            mp.loc['A'] = a_lst[i]
            mp.loc['C'] = c_lst[i]
            mp.loc['G'] = g_lst[i]
            mp.loc['T'] = t_lst[i]

            mw: int = motif_width_lst[i]
            sc: int = site_counts_lst[i]

            mp = norm_motif(mp, mw, alphabet)
            mp = apply_pseudocount_meme(mp, pseudocount, sc, mw, bgs, alphabet)

            motif: Motif = Motif(mp, mw, alphabet, motifID_lst[i], 
                                 motifName_lst[i])
            motif.setBg(bgs)

            motif_lst.append(motif)
        # end for

        return motif_lst

    finally:
        in_mtf.close()  # close the file anyway
    
    # end try

# end read_MEME_motif()


def process_motif_for_logodds(motif: Motif) -> Motif:
    """Computes the log-odds values from a probability matrix of a given
    PWM motif.

    During the computation of the log-odds matrix is also computed the 
    corresponding P-value matrix, using the dynamic programming 
    algorithm presented in Staden, 1994.

    Parameters
    ----------
    motif : Motif
        DNA motif 
        
    Returns
    -------
    Motif
        Input DNA motif with the log-odds matrix
    """

    # get the log-odds
    motif_log_odds: pd.DataFrame
    motif_log_odds = compute_log_odds(motif.getMotif_matrix(), motif.getWidth(),
                                      motif.getBg(), motif.getAlphabet())
    motif.setMotif_scoreMatrix(motif_log_odds)

    # scale the log-odds scores
    scaled_scores: np.ndarray
    min_val: int
    max_val: int
    scale: int
    offset: np.double

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
    pval_mat: np.array
    pval_mat = comp_pval_mat(motif)
    motif.setMotif_pval_matrix(pval_mat)

    motif.setMotif_scoreMatrix(scaled_scores.values)

    return motif

# end of process_motif_for_logodds()


def scale_pwm(motif_matrix: pd.DataFrame,
              alphabet: List[str],
              motif_width: int
) -> Tuple[np.ndarray, int, int, int, np.double]:
    """Scale the log-odds values of the motif scoring matrix.

    The values are scaled in the range [0, 1000]. The scaling improves
    computational speed while computing the score for each motif 
    occurrence candidate, and allows a constant time computation of 
    the corresponding P-value. 
        
    Parameters
    ----------
    motif_matrix : pd.DataFrame
        motif log-odds matrix
    alphabet: list
        DNA motif alphabet
    motif_width: int
        motif width

    Returns
    -------
    numpy.ndarray
        scaled motif scoring matrix
    int
        minimum value of the scaled scoring matrix
    int
        maximum value of the scaled scoring matrix
    int
        scaling factor
    numpy.double
        scaling offset
    """

    errmsg: str
    if not isinstance(motif_matrix, pd.DataFrame):
        errmsg = "\n\nERROR: The given motif matrix must be an instance of pandas.DataFrame"
        raise NoDataFrameException(errmsg)

    if motif_matrix.empty:
        errmsg = "\n\nERROR: The given motif matrix is empty"
        raise NotValidMotifMatrixException(errmsg)

    if not isinstance(alphabet, list):
        errmsg = "\n\nERROR: The alphabet given is not in a list"
        raise NotValidAlphabetException(errmsg)

    if not isListEqual(alphabet, DNA_ALPHABET):
        errmsg = "\n\nERROR: The alphabet given is not a valid DNA alphabet"
        raise NotValidAlphabetException(errmsg)

    assert motif_width > 0

    min_val: int
    max_val: int
    motif_matrix_sc: pd.DataFrame

    min_val = min(motif_matrix.min())
    max_val = max(motif_matrix.max())
    motif_matrix_sc = pd.DataFrame(index=list(motif_matrix.index), 
                                   columns=list(motif_matrix.columns), data=0)

    lower: int = min_val
    upper: int = max_val

    if lower == upper:  # all values are equal
        lower = np.double(upper - 1)

    offset: np.double
    scale_factor: int

    lower = np.floor(lower)
    offset = np.round(np.floor(lower))
    scale_factor = np.floor(RANGE / (upper - lower))

    # values will be in [0, 1000]
    for nuc in alphabet:
        for j in range(motif_width):
            scaled_score = np.round(
                (motif_matrix.loc[nuc, j] - (offset)) * scale_factor
                                   )
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


def get_motif_pwm(motif_file: str,
                  args_obj: Findmotif,
                  cores: int
) -> List[Motif]:
    """Starting point for the construction of a Motif object.

    The motif PWM will be read accordingly to the file format. From the
    read data will be computed the motif scoring matrix (with scores
    scaled) and the corresponding P-value matrix.

    All these data will be stored in a new Motif object.

    Parameters
    ----------
    motif_file : str
        path to the motif PWM
    args_obj : Findmotif
        container for arguments needed for the motif scoring and 
        P-value matrix computations
    cores : int
        number of cores to use during the computation (used only when
        processing MEME motif files)
    
    Returns
    -------
    List[Motif]
        processed Motif object as element of a list
    """

    bgs: dict
    pseudo: float
    no_reverse: bool
    verbose: bool

    # get arguments required to process the motif
    bgs = args_obj.get_bgfile()
    pseudo = args_obj.get_pseudo()
    no_reverse = args_obj.get_no_reverse()
    verbose = args_obj.get_verbose()

    errmsg: str
    if not motif_file:
        errmsg = "\n\nERROR: the motif file is missing"
        raise FileNotFoundError(errmsg)

    if (not isMEME_ff(motif_file)) and (not isJaspar_ff(motif_file)):
        errmsg = "\n\nERROR: the motif file must be in MEME or JASPAR format"
        raise NotValidFFException(errmsg)

    if isJaspar_ff(motif_file):
        motif = build_motif_JASPAR(motif_file, bgs, pseudo, no_reverse, verbose)

    elif isMEME_ff(motif_file):
        motif = build_motif_MEME(motif_file, bgs, pseudo, no_reverse, cores, 
                                 verbose)

    else:
        errmsg = ' '.join(["\n\nERROR: do not know what to do with file", 
                           motif_file])
        raise NotValidFFException(errmsg)
    # end if

    if not isinstance(motif, list):
        motif = [motif]

    return motif

# end of get_motif_pwm()


def pseudo_bg(bgs: Dict, 
              no_reverse: bool
) -> Dict:
    """Add a pseudocount and normalize the background probabilities of 
    nucleotides used to build the motif scoring matrix.

    A pseudocount value is added to the background probability 
    distribution.

    If are to be considered both the forward and the reverse strand the 
    background probabilities are averaged for the two strands. 

    The resulting background probabilities are then normalized.

    Parameters
    ----------
    bgs : dict
        background probability distribution
    no_reverse : bool
        if False only the forward strand will be considered, otherwise
        both forward and reverse are considered

    Returns 
    -------
    dict
        normalized background probablity distribution
    """

    bgs_avg: Dict
    bgs_proc: Dict

    errmsg: str
    if not isinstance(bgs, dict):
        errmsg = "\n\nERROR: unable to add the pseudocount to the background"
        raise NotValidBGException(errmsg)

    if not isinstance(no_reverse, bool):
        errmsg = ' '.join(["Boolean value required, got", str(type(no_reverse))])
        raise ValueException(errmsg)

    if not no_reverse:
        bgs_avg = average_bg_with_rc(bgs)
    else:
        bgs_avg = bgs

    bgs_proc = norm_bg(bgs_avg)

    return bgs_proc

# end pseudo_bg()


def average_bg_with_rc(bgs):
    """The background probabilities are averaged with the frequencies of
    the reverse complement strand.

    This step is executed only if the reverse complement has to be 
    considered in GRAFIMO analysis.

    Parameters
    ----------
    bgs : dict
        background probability distribution

    Returns
    -------
    dict
        background probability distribution averaged for reverse 
        complement feequencies 
    """

    bgs_avg = dict()

    for nuc in bgs.keys():
        rc: str = REV_COMPL[nuc]

        if REV_COMPL[rc] == nuc and ord(nuc) < ord(rc):
            avg_freq = np.double((bgs[nuc] + bgs[rc]) / np.double(2))
            bgs_avg.update({nuc: avg_freq})
            bgs_avg.update({rc: avg_freq})
        # end if
    # end for

    return bgs_avg

# end average_bg_with_rc()


def norm_bg(bgs):
    """Normalize the background probability distribution

    Parameters
    ----------
    bgs : dict
        background probability distribution
    
    Returns
    -------
    dict
        normalized background probability distribution
    """

    alphabet: List[str]
    tot: np.double

    # PSEUDO = np.double(0.0000005)

    alphabet = sorted(list(bgs.keys()))
    tot = np.double(len(alphabet) * PSEUDObg)
    bgs_norm = dict()

    for nuc in bgs.keys():
        tot += np.double(bgs[nuc])

    assert tot > 0

    prob: np.double
    for nuc in bgs.keys():
        prob = np.double((bgs[nuc] + PSEUDObg) / tot)
        bgs_norm.update({nuc: prob})

    tot = np.double(0)
    for nuc in bgs.keys():
        tot += bgs[nuc]

    assert tot != 0

    return bgs_norm

# end norm_bg()


def norm_motif(motif_probs: pd.DataFrame, 
               motif_width: int, 
               alphabet: List[str]
) -> pd.DataFrame:
    """Normalize motif PWM frequencies.

    The values normalized are given as probabilites rather than simple
    raw counts.

    Parameters
    ----------
    motif_probs : pandas.DataFrame
        motif probbaility matrix (PWM)
    motif_width : int
        motif width
    alphabet : list
        DNA motif alphabet

    Returns
    -------
    pandas.DataFrame
        normalized motif probability matrix
    """

    # allowed tolerance in the difference between the position 
    # probability and 1
    tolerance: float = 0.00001
    tot: np.double  

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

