"""Functions to process the motif Position Weight Matrix.

The PWM is processed to get log-odds score matrix (scaled in [0,1000]) and the 
corresponding P-value matrices, used during potential motif occurrences
scoring.

The resulting scoring and P-value matrices are stored in the corresponding, 
Motif objects.

TODO: manage TRANSFAC and PFM motif formats.

"""


from grafimo.GRAFIMOException import NoDataFrameException, ValueException, \
    NotValidMotifMatrixException, BGFileError, MotifFileReadError, \
    NotValidAlphabetException, NotValidFFException, MotifFileFormatError
from grafimo.utils import die, DNA_ALPHABET, REV_COMPL, PSEUDObg, isListEqual,\
     isJaspar_ff, isMEME_ff, almost_equal, RANGE, printProgressBar, \
     sigint_handler, exception_handler, UNIF
from motif_processing import readBGfile, get_uniformBG, \
    apply_pseudocount_jaspar, apply_pseudocount_meme, compute_log_odds, \
    comp_pval_mat
from grafimo.workflow import Findmotif
from grafimo.motif import Motif

from typing import List, Optional, Dict, Tuple

import multiprocessing as mp
import pandas as pd
import numpy as np

import signal
import time
import os


def build_motif_JASPAR(
    motif_file: str,
    bg_file: str,
    pseudocount: float,
    no_reverse: bool,
    verbose: bool, 
    debug: bool
) -> Motif:
    """Build the Motif object from a JASPAR motif Position Weight
    Matrix.

    It is computed the scoring matrix from the values given with the PWM
    and the P-value matrix to assign a statistical significance to
    each motif occurrence candidate, based on the resulting log-odds
    score.

    ...

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
    debug : bool
        trace the full error stack

    Returns
    -------
    Motif
        processed motif object
    """

    if not isinstance(motif_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_file).__name__), debug)
    if not os.path.isfile(motif_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(motif_file), debug)
    if not isJaspar_ff(motif_file, debug):
        errmsg = "Required JASPAR motif PWM parsing, but {} is not in JASPAR format.\n"
        exception_handler(MotifFileFormatError, errmsg.format(motif_file), debug)
    if not isinstance(bg_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(bg_file).__name__), debug)
    if bg_file != UNIF and not os.path.isfile(bg_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(bg_file), debug)
    if pseudocount <= 0:
        errmsg = "Pseudocount value must be positive.\n"
        exception_handler(ValueError, pseudocount, debug)
    if not isinstance(no_reverse, bool):
        errmsg = "Expected bool, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(no_reverse).__name__), debug)
    
    # parse motif PWM
    motif: Motif = read_JASPAR_motif(
        motif_file, bg_file, pseudocount, no_reverse, verbose, debug
    )
    if verbose: start_mp: float = time.time()
    motif = process_motif_for_logodds(motif, debug)  #  get log-odds values for motif
    if verbose:
        end_mp: float = time.time()
        print(
            "Motif %s processed in %.2fs" % (motif.motifID, (end_mp - start_mp))
        )

    return motif

# end of build_motif_JASPAR()


def read_JASPAR_motif(
    motif_file: str,
    bg_file: str,
    pseudocount: float,
    no_reverse: bool,
    verbose: bool, 
    debug: bool
) -> Motif:
    """Read a motif PWM in JASPAR format.

    The data read are then used to build the scoring matrix for the 
    motif, the P-value matrix, etc.

    ...

    Parameters
    ----------
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
    debug:
        trace the full error stack

    Returns
    -------
    Motif
        Motif object 
    """

    nucs: List[str] = list()
    counts: List[float] = list()
    if verbose:
        start_rm: float = time.time()
    try:
        ifstream = open(motif_file, mode="r")
        readlines = 0  # check for empty files
        # begin parsing
        header: str = str(ifstream.readline().strip()[1:])  
        if not header:  # empty file?
            errmsg = "{} seems to empty.\n"
            exception_handler(IOError, errmsg.format(motif_file), debug) 
        motifID, motifName = header.split('\t')[0:2]  
        readlines += 1
        while True:
            line = ifstream.readline().strip()
            if not line: break  # EOF or empty file?
            nuc = line.strip()[:1]
            count = list(map(float, line.strip()[1:].split()[1:][:-1]))
            nucs.append(nuc.upper())
            counts.append(count)
            readlines += 1
        if readlines <= 1:  # only header read ?
            errmsg = "{} seems to be empty.\n"
            exception_handler(IOError, errmsg.format(motif_file), debug)
    except:
        errmsg = "An error occurred while reading {}.\n"
        exception_handler(MotifFileReadError, errmsg.format(motif_file), debug)
    else:
        if any([len(c) != len(counts[0]) for c in counts]):
            errmsg = "Motif counts width mismatch.\n"
            exception_handler(ValueError, errmsg, debug)
        nucsmap = dict()  # used with np object
        for i in range(len(nucs)): nucsmap.update({nucs[i]:i})
        motif_counts: pd.DataFrame = pd.DataFrame(data=counts, index=nucs)  # motif count matrix
        motif_width: int = int(len(counts[0]))  
        alphabet: list = sorted(nucs) 
        
        # compute background
        if bg_file == UNIF: bgs = get_uniformBG(alphabet, debug)
        elif os.path.isfile(bg_file): bgs = readBGfile(bg_file, debug)
        else:
            errmsg = "Unable to parse {}.\n"
            exception_handler(BGFileError, errmsg.format(bg_file), debug)
        bgs = pseudo_bg(bgs, no_reverse, debug)  # add pseudocount to bg
        
        # motif probability matrix
        motif_probs = (motif_counts / motif_counts.sum(0)) 
        motif_probs = norm_motif(motif_probs, motif_width, alphabet, debug)
        motif_probs = apply_pseudocount_jaspar(
            motif_counts.to_numpy(), motif_probs.to_numpy(), pseudocount, bgs, 
            motif_width, alphabet, nucsmap, debug
        )
        motif: Motif = Motif(
            motif_probs, motif_width, alphabet, motifID, motifName, nucsmap
        )
        motif.setBg(bgs)

        if verbose:
            end_rm: float = time.time()
            msg: str = "Read motif %s in %.2fs" % (motifID, (end_rm - start_rm))
            print(msg)
    finally:
        ifstream.close() 

    return motif

# end of read_JASPAR_motif()


def build_motif_MEME(
    motif_file: str,
    bg_file: str,
    pseudocount: float,
    no_reverse: bool,
    cores: int,
    verbose: bool,
    debug: bool
) -> List[Motif]:
    """Read motif PWMs in MEME format.

    It is computed the scoring matrix from the values given with the PWM
    and the P-value matrix to assign a statistical significance to
    each motif occurrence candidate, based on the resulting log-odds
    score.

    ...

    Parameters:
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
    cores : int
        number of CPU cores (used when MEME file has more than one PWM)
    verbose : bool
        print additional information
    debug : bool
        trace the full error stack

    Returns
    -------
    Motif
        Motif object storing the data contained in motif_file
    """

    if not isinstance(motif_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_file).__name__), debug)
    if not os.path.isfile(motif_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(motif_file), debug)
    if not isMEME_ff(motif_file, debug):
        errmsg = "Required MEME motif PWM parsing, but {} is not in MEME format.\n"
        exception_handler(MotifFileFormatError, errmsg.format(motif_file), debug)

    if verbose: start_rm_all: float = time.time()
    motif_lst: List[Motif] = read_MEME_motif(
        motif_file, bg_file, pseudocount, no_reverse, verbose, debug
    )
    motif_num: int = len(motif_lst)
    if verbose:
        end_rm_all: float = time.time()
        print(
            "Read all motifs in %s in %.2fs." % (motif_file, (end_rm_all - start_rm_all))
        )
    print("\nRead {} motifs in {}".format(motif_num, motif_file))
    print("\nProcessing motifs\n")

    complete_motifs = list()  # fully processed motifs
    if verbose: start_mp_all: str = time.time()
    if motif_num >= cores:  # worth to use multiprocessing
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool: mp.Pool = mp.Pool(processes=cores)  
        # overwrite the default SIGINT handler to exit gracefully
        # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
        signal.signal(signal.SIGINT, original_sigint_handler)  

        try:
            args = [(motif, debug) for motif in motif_lst]
            res = (pool.starmap_async(process_motif_for_logodds, args))
            it: int = 0
            # ---- progress bar
            while (True):
                if res.ready():
                    # when finished call for the last time printProgressBar()
                    printProgressBar(
                        tot, tot, prefix='Progress:', suffix='Complete', length=50
                    )
                    break
                if it == 0: tot = res._number_left
                remaining = res._number_left
                printProgressBar(
                    (tot - remaining), tot, prefix='Progress:', suffix='Complete', length=50
                )
                time.sleep(1)
                it += 1
            complete_motifs += res.get(60 * 60 * 60)  # does not ignore signals
        except KeyboardInterrupt:
            pool.terminate()
            sigint_handler()
        else:
            pool.close()
            if verbose:
                end_mp_all: float = time.time()
                print(
                    "Processed motif(s) in %s in %.2fs" % (motif_file, (end_mp_all - start_mp_all))
                )
            return complete_motifs
    else: 
        for m in motif_lst:  # process each found motif
            complete_motifs.append(process_motif_for_logodds(m, debug))
        if verbose:
            end_mp_all: float = time.time()
            print(
                "Processed motif(s) in %s in %.2fs" % (motif_file, (end_mp_all - start_mp_all))
            )
        return complete_motifs

# end build_motif_MEME()


def read_MEME_motif(motif_file: str,
                    bg_file: str,
                    pseudocount: float,
                    no_reverse: bool,
                    verbose: bool,
                    debug: bool
) -> List[Motif]:
    """Read motif PWM in MEME format.

    The data read are then used to build the scoring matrix for the 
    motif, the P-value matrix, etc.

    Since a MEME file can contain one or more motifs, for each stored PWM
    is built the corresponding Motif object. The resulting set of motifs are 
    stored in a list, which will constitute a MotifSet object.

    ...
    
    Parameters
    ----------
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
    debug:
        trace the full error stack

    Returns
    -------
    List[Motif]
        list of Motif objects
    """

    if not isinstance(motif_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_file).__name__), debug)
    if not os.path.isfile(motif_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(motif_file), debug)
    if not isinstance(bg_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(bg_file).__name__), debug)
    if bg_file != UNIF and not os.path.isfile(bg_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(bg_file), debug)
    if not isinstance(pseudocount, float):
        errmsg = "Expected float, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(pseudocount).__name__), debug)
    if pseudocount <= 0:
        errmsg = "The pseudocount must be > 0.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(no_reverse, bool):
        errmsg = "Expected bool, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(no_reverse).__name__), debug)

    motifs_raw = list()
    motifs: List[Motif] = list()
    motifs_num = 0
    proceed = False
    # begin motif parsing
    try:
        ifstream = open(motif_file, mode="r")
        alphabet = __read_alphabet_meme(motif_file, ifstream, debug)  # shared by all motifs 
        nucsmap = dict()  # used with np object
        for i in range(len(alphabet)): nucsmap.update({alphabet[i]:i})
        while True:
            for line in ifstream:
                if line.startswith("MOTIF"): break  # new motif instance
            else:
                assert motifs_num == len(motifs_raw)
                proceed = True
                break
            if proceed: break  # read all motifs
            if verbose: start_rm = time.time()
            motifids = line.split()
            if len(motifids) == 2:  # only name
                motif_id = motifids[1]
                motif_name = motif_id
            else:  # assume first two fieds: id, name
                motif_id, motif_name = motifids[1:3]
            statistics = __read_statistics_meme(motif_file, ifstream, debug)
            probs = __read_counts_meme(
                motif_file, ifstream, statistics["width"], debug
            )
            motifs_raw.append(
                {"motifId":motif_id, "motifName":motif_name, "statistics":statistics, "counts":probs}
            )
            motifs_num += 1
            if verbose:
                end_rm = time.time()
                print("Read motif %s in %.2fs." % (motif_name, (end_rm - start_rm)))
        if not proceed:
            errmsg = "Unexpected premature EOF in {}.\n"
            exception_handler(EOFError, errmsg.format(motif_file), debug)
    except:
        errmsg = "An error occurred while reading {}.\n"
        exception_handler(MotifFileReadError, errmsg.format(motif_file), debug)
    else:
        if bg_file == UNIF: bgs = get_uniformBG(alphabet, debug)
        elif os.path.isfile(bg_file): bgs = readBGfile(bg_file, debug)
        else:
            errmsg = "Unable to parse {}.\n"
            exception_handler(BGFileError, errmsg.format(bg_file), debug)
        bgs = pseudo_bg(bgs, no_reverse, debug)  # add pseudocount to bg
        for i in range(motifs_num):
            mp = pd.DataFrame(np.matrix(motifs_raw[i]["counts"]))
            mp.index = alphabet
            mp = norm_motif(mp, motifs_raw[i]["statistics"]["width"], alphabet, debug)
            mp = apply_pseudocount_meme(
                mp.to_numpy(), pseudocount, motifs_raw[i]["statistics"]["nsites"],
                motifs_raw[i]["statistics"]["width"], bgs, alphabet, nucsmap, debug
            )
            motif: Motif = Motif(
                mp, motifs_raw[i]["statistics"]["width"], alphabet, 
                motifs_raw[i]["motifId"], motifs_raw[i]["motifName"], nucsmap
            )
            motif.setBg(bgs)
            motifs.append(motif)
    finally:
        ifstream.close()
    
    return motifs

# end read_MEME_motif()


def __read_alphabet_meme(motif_file: str, ifstream, debug: bool) -> List[str]:
    """Read alphabet from MEME files.
    
    ...

    Parameters
    ----------
    motif_file : str
        path to motif PWM
    ifstream : _io.TextIOWrapper
        input stream
    debug : bool
        trace the full error stack

    Returns
    -------
    list
        alphabet
    """
    
    for line in ifstream:
        if line.startswith("ALPHABET"): break
    else:
        errmsg = "Unexpected EOF reached, unable to parse {}.\n"
        exception_handler(EOFError, errmsg.format(motif_file), debug)
    if not line.startswith("ALPHABET"):
        errmsg = "No line stores alphabet in {}.\n"
        exception_handler(ValueError, errmsg.format(motif_file), debug)
    line = line.strip().replace("ALPHABET= ", "")
    if line == "ACGT": alphabet = sorted(list(line))
    else:
        errmsg = "The motif is not built on DNA alphabet.\n"
        exception_handler(ValueError, errmsg, debug)
    assert isListEqual(alphabet, DNA_ALPHABET)
    return alphabet


def __read_statistics_meme(motif_file: str, ifstream, debug: bool) -> Dict:
    """Read motif statistics form MEME files.

    ...

    Parameters
    ----------
    motif_file : str
        path to motif PWM
    ifstream : _io.TextIOWrapper
        input stream
    debug : bool
        trace the full error stack

    Returns 
    -------
    dict
        motif statistics
    """

    for line in ifstream:
        if line.startswith("letter-probability matrix:"): break  # statistics start here
    width = int(line.split("w=")[1].split()[0])
    nsites = int(line.split("nsites=")[1].split()[0])
    evalue = np.double(line.split("E=")[1].split()[0])
    statistics = {"width":width, "nsites":nsites, "evalue":evalue}
    return statistics


def __read_counts_meme(
    motif_file: str, 
    ifstream, 
    width: int, 
    debug: bool
) -> List[List[np.double]]:
    """Read motif letter probabilities from MEME files.

    ...

    Parameters
    ----------
    motif_file : str
        path to motif PWM
    ifstream : _io.TextIOWrapper
        input stream
    width : int
        motif width
    debug:
        trace the full error stack

    Returns 
    -------
    list
        motif letter probabilities
    """

    a = list()
    c = list()
    g = list()
    t = list()
    pos = 0
    for line in ifstream:
        freqs = line.split()
        if len(freqs) != 4:
            if pos < width:
                errmsg = "Unexpected end of motif found.\n"
                exception_handler(EOFError, errmsg, debug)
            break  # motif stop
        a.append(np.double(freqs[0]))
        c.append(np.double(freqs[1]))
        g.append(np.double(freqs[2]))
        t.append(np.double(freqs[3]))
        pos += 1
    probs = [a, c, g, t]
    if any([len(p) != len(probs[0]) for p in probs]):
        errmsg = "Mismatch in letter probabilities vectors lengths.\n"
        exception_handler(ValueError, errmsg, debug)
    return probs


def process_motif_for_logodds(motif: Motif, debug: bool) -> Motif:
    """Computes log-odds values from motif probability matrix (PFM).

    While processing  motif probability matrix for log-odds values is also
    computed the p-value matrix for the current motif PWM. 

    ...

    Parameters
    ----------
    motif : Motif
        DNA motif 
    debug : bool
        trace the full error stack
        
    Returns
    -------
    Motif
        motif log-odds matrix
    """

    if not isinstance(motif, Motif):
        errmsg = "Expected Motif, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif).__name__), debug)

    # compute log-odds
    motif_log_odds = compute_log_odds(
        motif.countMatrix, motif.width, motif.bg, motif.alphabet, motif.nucsmap, 
        debug
    )
    motif.set_motifScoreMatrix(motif_log_odds)

    # log-odds matrix scaling
    scaled_scores, min_val, max_val, scale, offset = scale_pwm(
        motif.scoreMatrix, motif.alphabet, motif.width, motif.nucsmap, debug
    )
    motif.set_motifScoreMatrix(scaled_scores)
    motif.set_isScaled()
    motif.set_scale(scale)
    motif.set_minVal(min_val)
    motif.set_maxVal(max_val)
    motif.set_offset(offset)

    # compute p-value matrix
    pval_mat = comp_pval_mat(motif, debug)
    motif.set_motifPvalMatrix(pval_mat)

    return motif

# end of process_motif_for_logodds()


def scale_pwm(
    motif_matrix: np.ndarray, 
    alphabet: List[str],
    motif_width: int,
    nucsmap: dict,
    debug: bool
) -> Tuple[np.ndarray, int, int, int, np.double]:
    """Scale the motif log-odds matrix scores to integer values.

    The values are scaled in the range [0, 1000]. The scaling improves
    computational speed while scoring potential motif occurrences, and allows
    constant time p-value estimatimation.

    ...
        
    Parameters
    ----------
    motif_matrix : numpy.ndarray
        motif log-odds matrix
    alphabet: list
        DNA motif alphabet
    motif_width: int
        motif width
    nucsmap: dict
        nucleotide index map
    debug : bool
        trace the full error stack

    Returns
    -------
    numpy.ndarray
        scaled motif score matrix
    int
        minimum value of the scaled score matrix
    int
        maximum value of the scaled score matrix
    int
        scaling factor
    numpy.double
        scaling offset
    """

    if not isinstance(motif_matrix, np.ndarray):
        errmsg = "Expected numpy.ndarray, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_matrix).__name__), debug)
    if motif_matrix.size == 0 or sum(sum(motif_matrix)) == 0:
        errmsg = "The motif log-odds natrix is empty.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(alphabet, list):
        errmsg = "Expected list, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(alphabet).__name__), debug)
    if not isListEqual(alphabet, DNA_ALPHABET):
        errmsg = "The motif is not built on DNA alphabet.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(motif_width, int):
        errmsg = "Expected int, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_width).__name__), debug)
    if motif_width <= 0:
        errmsg = "Forbidden motif width.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(nucsmap, dict):
        errmsg = "Expected dict, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(nucsmap).__name__), debug)

    min_val = motif_matrix.min()
    max_val = motif_matrix.max()
    motif_matrixsc = np.zeros(motif_matrix.shape, dtype=np.double)
    
    lower: int = min_val
    upper: int = max_val
    if lower == upper:  # all values are equal
        lower = np.double(upper - 1)
    lower = np.floor(lower)
    offset = np.round(np.floor(lower))
    scale_factor = np.floor(RANGE / (upper - lower))
    # values scaled in [0, 1000]
    for nuc in alphabet:
        for j in range(motif_width):
            scaled_score = np.round(
                (motif_matrix[nucsmap[nuc], j] - (offset)) * scale_factor
            )
            motif_matrixsc[nucsmap[nuc], j] = scaled_score
    # make sure the values are integers
    motif_matrixsc = motif_matrixsc.astype(int)
    min_val = int(motif_matrixsc.min())  # scaled min
    max_val = int(motif_matrixsc.max())  # scaled max

    return motif_matrixsc, min_val, max_val, int(scale_factor), offset

# end of scale_pwm()


def get_motif_pwm(motif_file: str,
                  args_obj: Findmotif,
                  cores: int,
                  debug: bool
) -> List[Motif]:
    """Construction of Motif object from PWM file.

    The motif PWM is processed in order to obtain the corresponding scoring
    matrix (values scaled in [0,1000]) and the corresponding P-value matrix,
    which is used to assign statistical significance to motif occurrence
    candidates scores.

    To store all these informations is created a Motif object.

    ...

    Parameters
    ----------
    motif_file : str
        path to motif PWM file (MEME or JASPAR format)
    args_obj : Findmotif
        arguments container
    cores : int
        CPU cores to use during motif processing (used only when
        processing MEME motif files with multiple PWMs)
    debug : bool
        trace the full error stack
    
    Returns
    -------
    List[Motif]
        Motif objects
    """

    bgs: dict = args_obj.bgfile
    pseudo: float = args_obj.pseudo
    no_reverse: bool = args_obj.noreverse
    verbose: bool = args_obj.verbose
    errmsg: str
    if not isinstance(motif_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_file).__name__), debug)
    if not os.path.isfile(motif_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(motif_file), debug)
    if (not isMEME_ff(motif_file, debug)) and (not isJaspar_ff(motif_file, debug)):
        errmsg = "Motif PWM must be in MEME or JASPAR format.\n"
        exception_handler(MotifFileFormatError, errmsg, debug)

    # chhose motif PWM parsing method
    if isJaspar_ff(motif_file, debug):
        motif = build_motif_JASPAR(
            motif_file, bgs, pseudo, no_reverse, verbose, debug
        )
    elif isMEME_ff(motif_file, debug):
        motif = build_motif_MEME(
            motif_file, bgs, pseudo, no_reverse, cores, verbose, debug
        )
    else:
        errmsg = "Motif PWM must be in MEME or JASPAR format.\n"
        exception_handler(MotifFileFormatError, errmsg, debug)

    # list instance required to proceed    
    if not isinstance(motif, list): motif = [motif]
    assert isinstance(motif, list) 
    return motif

# end of get_motif_pwm()


def pseudo_bg(bgs: Dict, 
              no_reverse: bool,
              debug: bool
) -> Dict:
    """Add pseudocount and normalize the nucleotides background probabilities. 
    The processed background probabilities are then used to compute the scoring
    matrix from the input motif PWM.

    When considered both forward and reverse strand, the background probabilities
    are weighted and averaged on both strands. 

    After the weighting and averaging steps (if required), the background 
    probabilities are normalized. 

    ...

    Parameters
    ----------
    bgs : dict
        background probability distribution
    no_reverse : bool
        if True the averaging and weighting operation on bg probabilities on 
        both DNA strands are skipped (only fwd strand considered). 
    debug: bool
        trace the full error stack

    Returns 
    -------
    dict
        normalized background probability distribution
    """

    if not isinstance(bgs, dict):
        errmsg = "Expected dict, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(bgs).__name__), debug)
    if not isinstance(no_reverse, bool):
        errmsg = "Expected bool, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(no_reverse).__name__), debug)
    
    if not no_reverse:  # fwd + rev strand
        bgs_avg = average_bg_with_rc(bgs, debug)
    else:  # only fwd
        bgs_avg = bgs
    bgs_proc = norm_bg(bgs_avg, debug)
    return bgs_proc

# end pseudo_bg()


def average_bg_with_rc(bgs: Dict, debug: bool):
    """Background probabilities are averaged with those occurring on the reverse
    complement strand.

    Parameters
    ----------
    bgs : dict
        background probability distribution
    debug: bool
        trace full error stack

    Returns
    -------
    dict
        background probability distribution averaged for reverse 
        complement feequencies 
    """

    if not isinstance(bgs, dict):
        errmsg = "Expected dict, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(bgs).__name__), debug)
    
    bgs_avg: Dict = dict()
    for nuc in bgs.keys():
        rc: str = REV_COMPL[nuc.upper()]
        if REV_COMPL[rc] == nuc and ord(nuc) < ord(rc):
            avg_freq = np.double((bgs[nuc] + bgs[rc]) / np.double(2))
            bgs_avg.update({nuc: avg_freq})
            bgs_avg.update({rc: avg_freq})
    return bgs_avg

# end average_bg_with_rc()


def norm_bg(bgs: Dict, debug: bool):
    """Normalize the background probability distribution.

    Parameters
    ----------
    bgs : dict
        background probability distribution
    debug: bool
        trace the full error stack
    
    Returns
    -------
    dict
        normalized background probability distribution
    """
    
    if not isinstance(bgs, dict):
        errmsg = "Expected dict, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(bgs).__name__), debug)

    alphabet: List[str] = sorted(list(bgs.keys()))
    tot = np.double(len(alphabet) * PSEUDObg)
    bgs_norm = dict()
    # PSEUDO = np.double(0.0000005)

    for nuc in bgs.keys(): tot += np.double(bgs[nuc])
    assert tot > 0
    for nuc in bgs.keys():
        prob = np.double((bgs[nuc] + PSEUDObg) / tot)
        bgs_norm.update({nuc: prob})
    tot = np.double(0)
    for nuc in bgs.keys(): tot += bgs[nuc]
    assert tot != 0
    return bgs_norm

# end norm_bg()


def norm_motif(motif_probs: pd.DataFrame, 
               motif_width: int, 
               alphabet: List[str],
               debug: bool
) -> pd.DataFrame:
    """Normalize motif PWM. The PWM values must be given as probability (so called
    PFM), rather than simple raw counts.

    Parameters
    ----------
    motif_probs : pandas.DataFrame
        motif probability matrix (PFM)
    motif_width : int
        motif width
    alphabet : list
        DNA motif alphabet
    debug: bool
        trace the full error stack

    Returns
    -------
    pandas.DataFrame
        normalized motif probability matrix (nPFM)
    """

    if not isinstance(motif_probs, pd.DataFrame):
        errmsg = "Expected pandas.DataFrame, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_probs).__name__), debug)
    if not isinstance(motif_width, int):
        errmsg = "Expected int, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_width).__name__), debug)
    if motif_width <= 0:
        errmsg = "Forbidden motif width.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(alphabet, list):
        errmsg = "Expected list, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(alphabet).__name__), debug)
    if any([nuc not in DNA_ALPHABET for nuc in alphabet]):
        errmsg = "The motif is not built on DNA alphabet.\n"
        exception_handler(ValueError, errmsg, debug)

    # tolerance in the difference between the position probability and 1
    tolerance: float = 0.00001 
    for j in range(motif_width):
        tot = np.double(0)
        for nuc in alphabet: tot += motif_probs.loc[nuc, j]
        assert tot != 0
        if not almost_equal(1, tot, tolerance):
            for nuc in alphabet:
                motif_probs.loc[nuc, j] = np.double(motif_probs.loc[nuc, j] / tot)

    return motif_probs

# end norm_motif()

