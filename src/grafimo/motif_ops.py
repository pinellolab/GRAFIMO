"""Functions to process the motif Position Weight Matrix.

The PWM is processed to get log-odds score matrix (scaled in [0,1000]) and the 
corresponding P-value matrices, used during potential motif occurrences
scoring.

The resulting scoring and P-value matrices are stored in the corresponding, 
Motif objects.

TODO: manage TRANSFAC and PFM motif formats.

"""


from grafimo.grafimo_errors import (
    BGFileError, 
    MotifFileReadError, 
    MotifFileFormatError
)
from grafimo.utils import (
    die, 
    isListEqual,
    isJaspar_ff,
    isMEME_ff,
    almost_equal,
    printProgressBar, 
    sigint_handler, 
    exception_handler, 
    DNA_ALPHABET, 
    REV_COMPL, 
    PSEUDOBG, 
    RANGE, 
    UNIF
)
from motif_processing import (
    read_bg_file, 
    get_uniform_bg, 
    apply_pseudocount_jaspar, 
    apply_pseudocount_meme, 
    compute_log_odds, 
    comp_pval_mat
)
from grafimo.workflow import Findmotif
from grafimo.motif import Motif

from typing import List, Optional, Dict, Tuple
from _io import TextIOWrapper

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
    """Build Motif object from a JASPAR motif Position Weight
    Matrix.

    It is computed the scoring matrix from the values given with the PWM
    and the P-value matrix to assign a statistical significance to
    each motif occurrence candidate, based on the resulting log-odds
    score.

    ...

    Parameters
    ----------
    motif_file : str
        Path to motif PWM file
    bg_file
        Path to the background file (Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html)).
    pseudocount : float
        Value to add to motif PWM counts
    no_reverse : bool
        Skip the reverse strand
    verbose : bool
        Print additional information
    debug : bool
        Trace the full error stack

    Returns
    -------
    Motif
        Processed motif 
    """

    if not isinstance(motif_file, str):
        errmsg = f"Expected {str.__name__}, got {type(motif_file).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not os.path.isfile(motif_file):
        errmsg = f"Unable to locate {motif_file}.\n"
        exception_handler(FileNotFoundError, errmsg, debug)
    if not isJaspar_ff(motif_file, debug):
        errmsg = f"Required JASPAR motif PWM parsing, but {motif_file} is not in JASPAR format.\n"
        exception_handler(MotifFileFormatError, errmsg, debug)
    if not isinstance(bg_file, str):
        errmsg = f"Expected {str.__name__}, got {type(bg_file).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if bg_file != UNIF and not os.path.isfile(bg_file):
        errmsg = f"Unable to locate {bg_file}.\n"
        exception_handler(FileNotFoundError, errmsg, debug)
    if pseudocount <= 0:
        errmsg = "Pseudocount value must be positive.\n"
        exception_handler(ValueError, pseudocount, debug)
    if not isinstance(no_reverse, bool):
        errmsg = f"Expected bool, got {no_reverse}.\n"
        exception_handler(TypeError, errmsg, debug)
    # parse motif PWM
    motif = read_JASPAR_motif(
        motif_file, bg_file, pseudocount, no_reverse, verbose, debug
    )
    if verbose: 
        start_mp = time.time()
    motif = process_motif_for_logodds(motif, debug)  #  get log-odds values for motif
    if verbose:
        end_mp = time.time()
        print(
            "Motif %s processed in %.2fs" % (motif.motif_id, (end_mp - start_mp))
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
    """Read motif PWMs in JASPAR format. The read data are used to 
    construct the corresponding Motif object.

    ...

    Parameters
    ----------
    motif_file : str
        Path to motif PWM file (JASPAR format)
    bg_file
        Path to the background file (Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html)).
    pseudocount : float
        Value to add to motif PWM counts
    no_reverse : bool
        Skip the reverse strand
    verbose : bool
        Print additional information
    debug:
        Trace the full error stack

    Returns
    -------
    Motif
        Motif object 
    """

    nucs = list()
    counts = list()
    if verbose:
        start_rm = time.time()
    try:
        handle = open(motif_file, mode="r")
        readlines = 0  # check for empty files
        # begin parsing
        header = str(handle.readline().strip()[1:])  
        if not header:  # empty file?
            errmsg = f"{motif_file} seems to empty.\n"
            exception_handler(IOError, errmsg, debug) 
        motif_id, motif_name = header.split("\t")[0:2]  
        readlines += 1
        while True:
            line = handle.readline().strip()
            if not line: 
                break  # EOF or empty file?
            nuc = line.strip()[:1]
            count = list(map(float, line.strip()[1:].split()[1:][:-1]))
            nucs.append(nuc.upper())
            counts.append(count)
            readlines += 1
        if readlines <= 1:  # only header read ?
            errmsg = f"{motif_file} seems to be empty.\n"
            exception_handler(IOError, errmsg, debug)
    except:
        errmsg = f"An error occurred while reading {motif_file}.\n"
        exception_handler(MotifFileReadError, errmsg, debug)
    else:
        if any([len(c) != len(counts[0]) for c in counts]):
            errmsg = "Motif counts width mismatch.\n"
            exception_handler(ValueError, errmsg, debug)
        nucsmap = dict()  # used with np object
        for i in range(len(nucs)): 
            nucsmap.update({nucs[i]:i})
        motif_counts = pd.DataFrame(data=counts, index=nucs)  # count matrix
        motif_width = int(len(counts[0]))  
        alphabet = sorted(nucs) 
        # compute background
        if bg_file == UNIF: 
            bgs = get_uniformBG(alphabet, debug)
        elif os.path.isfile(bg_file): 
            bgs = read_bg_file(bg_file, debug)
        else:
            errmsg = f"Unable to parse {bg_file}.\n"
            exception_handler(BGFileError, errmsg, debug)
        bgs = pseudo_bg(bgs, no_reverse, debug)  # add pseudocount to bg
        # motif probability matrix
        motif_probs = (motif_counts / motif_counts.sum(0)) 
        motif_probs = norm_motif(motif_probs, motif_width, alphabet, debug)
        motif_probs = apply_pseudocount_jaspar(
            motif_counts.to_numpy(), 
            motif_probs.to_numpy(), 
            pseudocount, 
            bgs, 
            motif_width, 
            alphabet, 
            nucsmap, 
            debug
        )
        motif = Motif(
            motif_probs, motif_width, alphabet, motif_id, motif_name, nucsmap
        )
        motif.set_bg(bgs)
        if verbose:
            end_rm = time.time()
            msg = "Read motif %s in %.2fs" % (motif_id, (end_rm - start_rm))
            print(msg)
    finally:
        handle.close() 
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
    """Read motif PWMs in MEME format. The read data are used to construct
    the corresponding motif object.

    The motif object stores the scoring matrix, used to assign a likelihood
    to each potential motif occurrence. Moreover, the P-value matrix is 
    used to assign a statistical significance value to each computed 
    score.

    ...

    Parameters:
    motif_file : str
        Path to motif PWM file
    bg_file
        Path to the background file (Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html))
    pseudocount : float
        Value to add to motif PWM counts
    no_reverse : bool
        Skip reverse complement
    cores : int
        Number of CPU cores (used when MEME file has more than one PWM)
    verbose : bool
        Print additional information
    debug : bool
        Trace the full error stack

    Returns
    -------
    Motif
        Motif object 
    """

    if not isinstance(motif_file, str):
        errmsg = f"Expected {str.__name__}, got {type(motif_file).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not os.path.isfile(motif_file):
        errmsg = f"Unable to locate {motif_file}.\n"
        exception_handler(FileNotFoundError, errmsg, debug)
    if not isMEME_ff(motif_file, debug):
        errmsg = f"Required MEME motif PWM parsing, but {motif_file} is not in MEME format.\n"
        exception_handler(MotifFileFormatError, errmsg, debug)

    if verbose: 
        start_rm_all = time.time()
    motif_lst = read_MEME_motif(
        motif_file, bg_file, pseudocount, no_reverse, verbose, debug
    )
    motif_num = len(motif_lst)
    if verbose:
        end_rm_all = time.time()
        print(
            "Read all motifs in %s in %.2fs." % (
                motif_file, (end_rm_all - start_rm_all)
            )
        )
    print(f"\nRead {motif_num} motifs in {motif_file}")
    print("\nProcessing motifs\n")
    complete_motifs = list()  # fully processed motifs
    if verbose: 
        start_mp_all = time.time()
    if motif_num >= cores:  # worth to use multiprocessing
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool = mp.Pool(processes=cores)  
        # overwrite the default SIGINT handler to exit gracefully
        # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
        signal.signal(signal.SIGINT, original_sigint_handler)  
        try:
            args = [(motif, debug) for motif in motif_lst]
            res = (pool.starmap_async(process_motif_for_logodds, args))
            it = 0
            # ---- progress bar
            while (True):
                if res.ready():
                    # when finished call for the last time printProgressBar()
                    printProgressBar(
                        tot, 
                        tot, 
                        prefix="Progress:", suffix="Complete", length=50
                    )
                    break
                if it == 0: 
                    tot = res._number_left
                remaining = res._number_left
                printProgressBar(
                    (tot - remaining), 
                    tot, 
                    prefix="Progress:", 
                    suffix="Complete", 
                    length=50
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
                end_mp_all = time.time()
                print(
                    "Processed motif(s) in %s in %.2fs" % (
                        motif_file, (end_mp_all - start_mp_all)
                    )
                )
            return complete_motifs
    else: 
        for m in motif_lst:  # process each found motif
            complete_motifs.append(process_motif_for_logodds(m, debug))
        if verbose:
            end_mp_all = time.time()
            print(
                "Processed motif(s) in %s in %.2fs" % (
                    motif_file, (end_mp_all - start_mp_all)
                )
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
    """Read motif PWM in MEME file format.

    The data read are then used to build the scoring matrix for the 
    motif, the P-value matrix, etc.

    Since a MEME file can contain one or more motifs, for each stored PWM
    is built the corresponding Motif object. The resulting set of motifs are 
    stored in a list, which will constitute a MotifSet object.

    ...
    
    Parameters
    ----------
    motif_file : str
        Path to motif PWM file (MEME format)
    bg_file
        Path to the background file (Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html)).
    pseudocount : float
        Value to add to motif PWM counts
    no_reverse : bool
        Skip reverse strand
    verbose : bool
        Print additional information
    debug:
        Trace the full error stack

    Returns
    -------
    List[Motif]
        List of Motif objects
    """

    if not isinstance(motif_file, str):
        errmsg = f"Expected {str.__name__}, got {type(motif_file).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not os.path.isfile(motif_file):
        errmsg = f"Unable to locate {motif_file}.\n"
        exception_handler(FileNotFoundError, errmsg, debug)
    if not isinstance(bg_file, str):
        errmsg = f"Expected {str.__name__}, got {type(bg_file).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if bg_file != UNIF and not os.path.isfile(bg_file):
        errmsg = f"Unable to locate {bg_file}.\n"
        exception_handler(FileNotFoundError, errmsg, debug)
    if not isinstance(pseudocount, float):
        errmsg = f"Expected {float.__name__}, got {type(pseudocount).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if pseudocount <= 0:
        errmsg = "The pseudocount must be > 0.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(no_reverse, bool):
        errmsg = f"Expected {bool.__name__}, got {type(no_reverse).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    motifs_raw = list()
    motifs = list()
    motifs_num = 0
    proceed = False
    # begin motif parsing
    try:
        handle = open(motif_file, mode="r")
        alphabet = __read_alphabet_meme(motif_file, handle, debug)  # shared by all motifs 
        nucsmap = dict()  # used with np object
        for i in range(len(alphabet)): nucsmap.update({alphabet[i]:i})
        while True:
            for line in handle:
                if line.startswith("MOTIF"): 
                    break  # new motif instance
            else:
                assert motifs_num == len(motifs_raw)
                proceed = True
                break
            if proceed: 
                break  # read all motifs
            if verbose: 
                start_rm = time.time()
            motifids = line.split()
            if len(motifids) == 2:  # only name
                motif_id = motifids[1]
                motif_name = motif_id
            else:  # assume first two fieds: id, name
                motif_id, motif_name = motifids[1:3]
            statistics = __read_statistics_meme(motif_file, handle, debug)
            probs = __read_counts_meme(
                motif_file, handle, statistics["width"], debug
            )
            motifs_raw.append(
                {
                    "motifId":motif_id, 
                    "motifName":motif_name, 
                    "statistics":statistics, 
                    "counts":probs
                }
            )
            motifs_num += 1
            if verbose:
                end_rm = time.time()
                print(
                    "Read motif %s in %.2fs." % (motif_name, (end_rm - start_rm))
                )
        if not proceed:
            errmsg = f"Unexpected premature EOF in {motif_file}.\n"
            exception_handler(EOFError, errmsg, debug)
    except:
        errmsg = f"An error occurred while reading {motif_file}.\n"
        exception_handler(MotifFileReadError, errmsg, debug)
    else:
        if bg_file == UNIF: 
            bgs = get_uniform_bg(alphabet, debug)
        elif os.path.isfile(bg_file): 
            bgs = read_bg_file(bg_file, debug)
        else:
            errmsg = f"Unable to parse {bg_file}.\n"
            exception_handler(BGFileError, errmsg, debug)
        bgs = pseudo_bg(bgs, no_reverse, debug)  # add pseudocount to bg
        for i in range(motifs_num):
            mp = pd.DataFrame(np.matrix(motifs_raw[i]["counts"]))
            mp.index = alphabet
            mp = norm_motif(
                mp, motifs_raw[i]["statistics"]["width"], alphabet, debug
            )
            mp = apply_pseudocount_meme(
                mp.to_numpy(), 
                pseudocount, 
                motifs_raw[i]["statistics"]["nsites"],
                motifs_raw[i]["statistics"]["width"], 
                bgs, 
                alphabet, 
                nucsmap, 
                debug
            )
            motif = Motif(
                mp, 
                motifs_raw[i]["statistics"]["width"], 
                alphabet, 
                motifs_raw[i]["motifId"], 
                motifs_raw[i]["motifName"], 
                nucsmap
            )
            motif.set_bg(bgs)
            motifs.append(motif)
    finally:
        handle.close()
    return motifs

# end read_MEME_motif()


def __read_alphabet_meme(
    motif_file: str, handle: TextIOWrapper, debug: bool
) -> List[str]:
    """Read alphabet from MEME files.
    
    ...

    Parameters
    ----------
    motif_file : str
        Path to motif PWM file
    handle : TextIOWrapper
        Input stream
    debug : bool
        Trace the full error stack

    Returns
    -------
    List[str]
        Motif alphabet
    """
    
    for line in handle:
        if line.startswith("ALPHABET"): 
            break
    else:
        errmsg = f"Unexpected EOF reached, unable to parse {motif_file}.\n"
        exception_handler(EOFError, errmsg, debug)
    if not line.startswith("ALPHABET"):
        errmsg = f"No line stores alphabet in {motif_file}.\n"
        exception_handler(ValueError, errmsg, debug)
    line = line.strip().replace("ALPHABET= ", "")
    if line == "ACGT": 
        alphabet = sorted(list(line))
    else:
        errmsg = "The motif is not built on DNA alphabet.\n"
        exception_handler(ValueError, errmsg, debug)
    assert isListEqual(alphabet, DNA_ALPHABET)
    return alphabet


def __read_statistics_meme(
    motif_file: str, handle: TextIOWrapper, debug: bool
) -> Dict[str, int]:
    """Read motif statistics form MEME files.

    ...

    Parameters
    ----------
    motif_file : str
        Path to motif PWM file
    handle : TextIOWrapper
        Input stream
    debug : bool
        Trace the full error stack

    Returns 
    -------
    Dict[str, int]
        Motif statistics
    """

    for line in handle:
        if line.startswith("letter-probability matrix:"): 
            break  # statistics start here
    width = int(line.split("w=")[1].split()[0])
    nsites = int(line.split("nsites=")[1].split()[0])
    evalue = np.double(line.split("E=")[1].split()[0])
    statistics = {"width":width, "nsites":nsites, "evalue":evalue}
    return statistics


def __read_counts_meme(
    motif_file: str, 
    handle: TextIOWrapper, 
    width: int, 
    debug: bool
) -> List[List[np.double]]:
    """Read motif letter probabilities from MEME files.

    ...

    Parameters
    ----------
    motif_file : str
        Path to motif PWM file
    handle : TextIOWrapper
        Input stream
    width : int
        Motif width
    debug:
        Trace the full error stack

    Returns 
    -------
    List[List[np.double]]
        Motif letter probabilities
    """

    a = list()
    c = list()
    g = list()
    t = list()
    pos = 0
    for line in handle:
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
    """Computes log-odds from motif probability matrix (PFM).

    While processing  motif probability matrix for log-odds values is 
    also computed the p-value matrix for the current motif PWM. 

    ...

    Parameters
    ----------
    motif : Motif
        Motif object 
    debug : bool
        Trace the full error stack
        
    Returns
    -------
    Motif
        Motif log-odds matrix
    """

    if not isinstance(motif, Motif):
        errmsg = f"Expected {type(Motif).__name__}, got {type(motif).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    # compute log-odds
    motif_log_odds = compute_log_odds(
        motif.count_matrix, 
        motif.width, 
        motif.bg, 
        motif.alphabet, 
        motif.nucsmap, 
        debug
    )
    motif.set_motif_score_matrix(motif_log_odds)
    # log-odds matrix scaling
    scaled_scores, min_val, max_val, scale, offset = scale_pwm(
        motif.score_matrix, 
        motif.alphabet, 
        motif.width, 
        motif.nucsmap, 
        debug
    )
    motif.set_motif_score_matrix(scaled_scores)
    motif.set_is_scaled()
    motif.set_scale(scale)
    motif.set_min_val(min_val)
    motif.set_max_val(max_val)
    motif.set_offset(offset)
    # compute p-value matrix
    pval_mat = comp_pval_mat(motif, debug)
    motif.set_motif_pval_matrix(pval_mat)
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
    computational speed while scoring potential motif occurrences, and 
    allows O(1) p-value estimatimation.

    ...
        
    Parameters
    ----------
    motif_matrix : numpy.ndarray
        Motif log-odds matrix
    alphabet: list
        Motif alphabet
    motif_width: int
        Motif width
    nucsmap: dict
        Nucleotide index map
    debug : bool
        Trace the full error stack

    Returns
    -------
    numpy.ndarray
        Scaled motif score matrix
    int
        Minimum value of the scaled score matrix
    int
        Maximum value of the scaled score matrix
    int
        Scaling factor
    numpy.double
        Scaling offset
    """

    if not isinstance(motif_matrix, np.ndarray):
        errmsg = f"Expected {type(np.ndarray).__name__}, got {type(motif_matrix).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if motif_matrix.size == 0 or sum(sum(motif_matrix)) == 0:
        errmsg = "The motif log-odds natrix is empty.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(alphabet, list):
        errmsg = f"Expected {list.__name__}, got {type(alphabet).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not isListEqual(alphabet, DNA_ALPHABET):
        errmsg = "The motif is not built on DNA alphabet.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(motif_width, int):
        errmsg = f"Expected {int.__name__}, got {type(motif_width).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if motif_width <= 0:
        errmsg = "Forbidden motif width.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(nucsmap, dict):
        errmsg = f"Expected {dict.__name__}, got {type(nucsmap).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    min_val = motif_matrix.min()
    max_val = motif_matrix.max()
    motif_matrixsc = np.zeros(motif_matrix.shape, dtype=np.double)
    lower = min_val
    upper = max_val
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


def get_motif_pwm(
    motif_file: str, args_obj: Findmotif, cores: int, debug: bool
) -> List[Motif]:
    """Construction of Motif object from PWM file.

    The motif PWM is processed in order to obtain the corresponding 
    scoring matrix (values scaled in [0,1000]) and the corresponding 
    P-value matrix, which is used to assign statistical significance to 
    motif occurrence candidates scores.

    To store all these informations is created a Motif object.

    ...

    Parameters
    ----------
    motif_file : str
        Path to motif PWM file (MEME or JASPAR format)
    args_obj : Findmotif
        Command line arguments container
    cores : int
        CPU cores to use during motif processing (used only when
        processing MEME motif files with multiple PWMs)
    debug : bool
        Trace the full error stack
    
    Returns
    -------
    List[Motif]
        Motif objects
    """

    bgs = args_obj.bgfile
    pseudo = args_obj.pseudo
    no_reverse = args_obj.noreverse
    verbose = args_obj.verbose
    if not isinstance(motif_file, str):
        errmsg = f"Expected {str.__name__}, got {type(motif_file).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not os.path.isfile(motif_file):
        errmsg = f"Unable to locate {motif_file}.\n"
        exception_handler(FileNotFoundError, errmsg, debug)
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
    if not isinstance(motif, list): 
        motif = [motif]
    assert isinstance(motif, list) 
    return motif

# end of get_motif_pwm()


def pseudo_bg(bgs: Dict, no_reverse: bool, debug: bool) -> Dict[str, float]:
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
    bgs : Dict[str, float]
        Background probability distribution
    no_reverse : bool
        Skip reverse strand 
    debug: bool
        Trace the full error stack

    Returns 
    -------
    Dict[str, float]
        Normalized background probability distribution
    """

    if not isinstance(bgs, dict):
        errmsg = f"Expected {dict.__name__}, got {type(bgs).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(no_reverse, bool):
        errmsg = f"Expected {bool.__name__}, got {type(no_reverse).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not no_reverse:  # fwd + rev strand
        bgs_avg = average_bg_with_rc(bgs, debug)
    else:  # only fwd
        bgs_avg = bgs
    bgs_proc = norm_bg(bgs_avg, debug)
    return bgs_proc

# end pseudo_bg()


def average_bg_with_rc(bgs: Dict, debug: bool) -> Dict[str, float]:
    """Background probabilities are averaged with those occurring on the 
    reverse strand.

    ...

    Parameters
    ----------
    bgs : Dict[str, float]
        Background probability distribution
    debug: bool
        Trace full error stack

    Returns
    -------
    Dict[str, float]
        Background probability distribution averaged for reverse 
        complement feequencies 
    """

    if not isinstance(bgs, dict):
        errmsg = f"Expected {dict.__name__}, got {type(bgs).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    bgs_avg = dict()
    for nuc in bgs.keys():
        rc = REV_COMPL[nuc.upper()]
        if REV_COMPL[rc] == nuc and ord(nuc) < ord(rc):
            avg_freq = np.double((bgs[nuc] + bgs[rc]) / np.double(2))
            bgs_avg.update({nuc: avg_freq})
            bgs_avg.update({rc: avg_freq})
    return bgs_avg

# end average_bg_with_rc()


def norm_bg(bgs: Dict, debug: bool) -> Dict[str, float]:
    """Normalize the background probability distribution.

    ...
    
    Parameters
    ----------
    bgs : Dict[str, float]
        Background probability distribution
    debug: bool
        Trace the full error stack
    
    Returns
    -------
    Dict[str, float]
        Normalized background probability distribution
    """
    
    if not isinstance(bgs, dict):
        errmsg = f"Expected {dict.__name__}, got {type(bgs).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    alphabet = sorted(list(bgs.keys()))
    tot = np.double(len(alphabet) * PSEUDOBG)
    bgs_norm = dict()
    # PSEUDO = np.double(0.0000005)
    for nuc in bgs.keys(): 
        tot += np.double(bgs[nuc])
    assert tot > 0
    for nuc in bgs.keys():
        prob = np.double((bgs[nuc] + PSEUDOBG) / tot)
        bgs_norm.update({nuc: prob})
    tot = np.double(0)
    for nuc in bgs.keys(): tot += bgs[nuc]
    assert tot != 0
    return bgs_norm

# end norm_bg()


def norm_motif(
    motif_probs: pd.DataFrame, 
    motif_width: int, 
    alphabet: List[str],
    debug: bool
) -> pd.DataFrame:
    """Normalize motif PWM. The PWM values must be given as probability 
    (so called PFM), rather than simple raw counts.

    ...

    Parameters
    ----------
    motif_probs : pandas.DataFrame
        Motif probability matrix (PFM)
    motif_width : int
        Motif width
    alphabet : list
        Motif alphabet
    debug: bool
        Trace the full error stack

    Returns
    -------
    pandas.DataFrame
        Normalized motif probability matrix (nPFM)
    """

    if not isinstance(motif_probs, pd.DataFrame):
        errmsg = f"Expected {pd.DataFrame}, got {type(motif_probs).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not isinstance(motif_width, int):
        errmsg = f"Expected {int.__name__}, got {type(motif_width).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if motif_width <= 0:
        errmsg = "Forbidden motif width.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(alphabet, list):
        errmsg = f"Expected {list.__name__}, got {type(alphabet).__name__}.\n"
        exception_handler(TypeError, errmsg, debug)
    if any([nuc not in DNA_ALPHABET for nuc in alphabet]):
        errmsg = "The motif is not built on DNA alphabet.\n"
        exception_handler(ValueError, errmsg, debug)
    # tolerance in the difference between the position probability and 1
    tolerance = 0.00001 
    for j in range(motif_width):
        tot = np.double(0)
        for nuc in alphabet: 
            tot += motif_probs.loc[nuc, j]
        assert tot != 0
        if not almost_equal(1, tot, tolerance):
            for nuc in alphabet:
                motif_probs.loc[nuc, j] = np.double(
                    motif_probs.loc[nuc, j] / tot
                )
    return motif_probs

# end norm_motif()

