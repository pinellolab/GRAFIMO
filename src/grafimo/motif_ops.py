"""Functions to process the motif Position Weight Matrix.

The PWM is processed to get log-odds score matrix (scaled in [0,1000]) and the 
corresponding P-value matrices, used during potential motif occurrences
scoring.

The resulting scoring and P-value matrices are stored in the corresponding, 
Motif objects.

TODO: manage TRANSFAC and PFM motif formats.

"""


from grafimo.GRAFIMOException import (
    BGFileError, 
    MotifFileReadError, 
    MotifFileFormatError
)
from grafimo.utils import (
    die, 
    isListEqual,
    isJaspar_ff, 
    isMEME_ff,
    isPFM_ff, 
    isTRANSFAC_ff,
    is_number,
    almost_equal,
    sigint_handler, 
    exception_handler,
    printProgressBar,
    DNA_ALPHABET, 
    REV_COMPL, 
    PSEUDObg, 
    RANGE,
    UNIF,
    JASPAR,
    MEME,
    PFM,
    TRANSFAC
)
from motif_processing import (
    readBGfile, 
    get_uniformBG, 
    apply_pseudocount_jaspar_pfm_transfac, 
    apply_pseudocount_meme, 
    compute_log_odds, 
    comp_pval_mat
)
from grafimo.workflow import Findmotif
from grafimo.motif import Motif

from typing import (
    List, 
    Optional, 
    Dict, 
    Tuple
)

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
    """Build the Motif object from JASPAR-like PFMs.

    Motif examples:

    # JASPAR
    >MA0001.1 AGL3
    A  [ 0  3 79 40 66 48 65 11 65  0 ]
    C  [94 75  4  3  1  2  5  2  3  3 ]
    G  [ 1  0  3  4  1  0  5  3 28 88 ]
    T  [ 2 19 11 50 29 47 22 81  1  6 ]
    
    or

    >MA0001.1 AGL3
    0  3 79 40 66 48 65 11 65  0
    94 75  4  3  1  2  5  2  3  3
    1  0  3  4  1  0  5  3 28 88
    2 19 11 50 29 47 22 81  1  6

    # scertf pcm
    A | 9 1 1 97 1 94
    T | 80 1 97 1 1 2
    C | 9 97 1 1 1 2
    G | 2 1 1 1 97 2


    Compute the scoring and the P-value matrix from PFM's values. The 
    P-value matrix is used to assign the statistical significance of 
    each binding affinity score obtained from the scoring matrix. The
    score of each potential motif occurremce is assigned as log-likelihood. 

    ...

    Parameters
    ----------
    motif_file : str
        path to motif PFM
    bg_file
        path to the background file in Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html).
    pseudocount : float
        value to add to PFM counts
    no_reverse : bool
        if False consider only the forward strand, forward and reverse
        otherwise
    verbose : bool
        print additional information
    debug : bool
        trace the full error stack

    Returns
    -------
    Motif
        Motif object
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
    if pseudocount <= 0:
        errmsg = "Pseudocount value must be positive.\n"
        exception_handler(ValueError, pseudocount, debug)
    if not isinstance(no_reverse, bool):
        errmsg = "Expected bool, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(no_reverse).__name__), debug)
    
    # parse motif PWM
    motif: Motif = read_JASPAR_motif(
        motif_file, 
        bg_file, 
        pseudocount, 
        no_reverse, 
        verbose, 
        debug
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
    """Parser for JASPAR-like motif files.

    ...

    Parameters
    ----------
    motif_file : str
        path to motif PFM
    bg_file
        path to the background file in Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html).
    pseudocount : float
        value to add to PFM counts
    no_reverse : bool
        if False consider only the forward strand, forward and reverse
        otherwise
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
    motif_id: str = ""
    motif_name: str = ""
    pwm_row_num = 0
    if verbose: start_rm = time.time()
    try:
        ifstream = open(motif_file, mode="r")
        while True:
            line = ifstream.readline().strip()
            if not line: break  # EOF or empty file?
            if line.startswith("#"): continue  # comment
            elif line.startswith(">"):
                header = line.split()
                if len(header) == 1:
                    motif_id = header[0]
                    motif_name = motif_id
                elif len(header) > 1:
                    motif_id = header[0]
                    motif_name = header[1]
                else:  # len(header) == 0  (????)
                    motif_id = motif_file.split("/")[-1]  # motif fname
                    motif_name = motif_id
            else:  # motif start
                columns = line.split()
                if columns[0] in DNA_ALPHABET:
                    nucs.append(columns[0])
                    columns = columns[1:]
                    if not is_number(columns[0], debug):
                        columns = columns[1:]
                        if not is_number(columns[-1], debug):
                            columns = columns[:-1]
                count = list(map(float, columns))
                counts.append(count)
                pwm_row_num += 1
            print(pwm_row_num)
            print(counts)
        if pwm_row_num <= 1:  # read only header or empty file
            errmsg = "{} seems to be empty.\n"
            exception_handler(IOError, errmsg.format(motif_file), debug)
        print(counts)
        assert pwm_row_num == 4
        if not motif_id: motif_id = motif_file.split("/")[-1]
        if not motif_name: motif_name = motif_file.split("/")[-1]
    except:
        errmsg = "An error occured while reading {}.\n"
        exception_handler(MotifFileReadError, errmsg.format(motif_file), debug)
    finally:
        ifstream.close()
    if any([len(c) != len(counts[0]) for c in counts]):
        errmsg = "Motif counts width mismatch.\n"
        exception_handler(ValueError, errmsg, debug)
    if not nucs: nucs = DNA_ALPHABET
    nucsmap = dict()  # used with np objects
    for i in range(len(nucs)): nucsmap[nucs[i]] = i
    # motif counts matrix
    motif_counts: pd.DataFrame = pd.DataFrame(data=counts, index=nucs)
    motif_width: int = motif_counts.shape[1]
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
    motif_probs = apply_pseudocount_jaspar_pfm_transfac(
        motif_counts.to_numpy(), 
        motif_probs.to_numpy(), 
        pseudocount, 
        bgs, 
        motif_width, 
        alphabet, 
        nucsmap, 
        debug
    )
    motif: Motif = Motif(
        motif_probs, 
        motif_width, 
        alphabet, 
        motif_id,
        motif_name,
        nucsmap
    )
    motif.setBg(bgs)
    if verbose:
        end_rm: float = time.time()
        msg: str = "Read motif %s in %.2fs" % (motif_id, (end_rm - start_rm))
        print(msg)

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
    if not isinstance(bg_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(bg_file).__name__), debug)
    if bg_file != UNIF and not os.path.isfile(bg_file):
        errmsg = "Unable to loacet {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(bg_file), debug)
    if pseudocount <= 0:
        errmsg = "Pseudocount must be > 0.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(no_reverse, bool):
        errmsg = "Expected boll, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(no_reverse).__name__), debug)

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
    if verbose: start_mp_all: float = time.time()
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

# end of __read_alphabet_meme()


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

# end of __read_statistics_meme()


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

# end of __read_counts_meme()


def build_PFM_motif(
    motif_file: str,
    bg_file: str,
    pseudocount: float,
    no_reverse: bool,
    cores: int,
    verbose: bool,
    debug: bool
) -> Motif:
    """Compute the motif score and P-value matrices from PFM motif files.

    File examples:

    # hocomoco
    >GATA1_HUMAN.H11MO.0.A
    96	202	103	99
    58	72	95	275
    28	39	320	113
    59	103	195	143
    97	157	141	105
    134	123	139	104
    112	105	184	99
    131	84	170	115
    107	112	157	124
    166	94	162	78
    92	194	161	53
    360	10	4	126
    3	0	497	0
    498	0	0	2
    1	0	6	493
    473	2	7	18
    429	8	45	18
    107	77	277	39
    176	94	173	57

    # tiffin
    T   A   G   C
    30  0   28  40
    0   0   0   99
    0   55  14  29
    0   99  0   0
    20  78  0   0
    0   52  7   39
    19  46  11  22
    0   60  38  0
    0   33  0   66
    73  0   25  0
    99  0   0   0

    Computes PSSM scoring matrix from counts values stored in the input 
    PFM file. While computing the PSSM, the function build the 
    corresponding P-value matrix. The PSSM assigns log-likelihood scores
    to potential motif occurrences, while the P-value matrix computes 
    scores statistical significance. 

    ...

    Parameters
    ----------
    motif_file : str
        path to motif file (PFM format)
    bg_file
        path to background file in Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html).
    pseudocount : float
        value to add to motif PWM counts
    no_reverse : bool
        if False, consider only forward strand. Use both strands 
        otherwise
    cores : int
        number of CPU cores 
    verbose : bool
        print additional information
    debug : bool
        trace the full error stack

    Returns
    -------
    Motif
        Motif object
    """
    
    if not isinstance(motif_file, str):
        errmsg = "\n\nExpected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_file).__name__), debug)
    if not os.path.isfile(motif_file):
        errmsg = "\n\nUnable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(motif_file), debug)
    if not isPFM_ff(motif_file, debug):
        errmsg = "Required PFM motif PWM parsing, but {} is not in PFM format.\n"
        exception_handler(MotifFileFormatError, errmsg.format(motif_file), debug)
    if not isinstance(bg_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(bg_file).__name__), debug)
    if bg_file != UNIF and not os.path.isfile(bg_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(bg_file), debug)
    if pseudocount <= 0:
        errmsg = "Pseudocount must be > 0.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(no_reverse, bool):
        errmsg = "Expected bool, got {}\n."
        exception_handler(TypeError, errmsg.format(type(no_reverse).__name__), debug)

    # parse motif
    if verbose: start_rm_all: float = time.time()
    motifs: List[Motif] = read_PFM_motif(
        motif_file, 
        bg_file, 
        pseudocount, 
        no_reverse, 
        verbose, 
        debug 
    )
    motif_num = len(motifs)
    if verbose:
        stop_rm_all = time.time()
        print(
            "Read all motifs in %s in %.2fs." % (motif_file, (stop_rm_all - start_rm_all))
        )
    print(f"\nRead {motif_num} motifs in {motif_file}")
    print("\nProcessing motifs\n")
    complete_motifs = list()
    if verbose: start_mp_all = time.time()
    if motif_num >= cores:  # worth using multiprocessing
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool: mp.Pool = mp.Pool(processes=cores)
        # overwrite the default SIGINT handler to exit gracefully
        # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
        signal.signal(signal.SIGINT, original_sigint_handler)
        try:
            args = [(motif, debug) for motif in motifs]
            res = (pool.starmap_async(process_motif_for_logodds, args))
            it: int = 0
            # --- progress bar
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
        pool.close()
        if verbose:
            stop_mp_all: float = time.time()
            print(
                "Processed motif(s) in %s in %.2fs" % (motif_file, (stop_mp_all - start_mp_all))
            )
        return complete_motifs
    else:
        for motif in motifs:
            complete_motifs.append(process_motif_for_logodds(motif, debug))
        if verbose:
            stop_mp_all: float = time.time()
            print(
                "Processed motif(s) in %s in %.2fs" % (motif_file, (stop_mp_all - start_mp_all))
            )
        return complete_motifs
    
# end of build_PFM_motif()


def read_PFM_motif(
    motif_file: str,
    bg_file: str,
    pseudocount: float,
    no_reverse: bool,
    verbose: bool,
    debug: bool
) -> Motif:
    """Read motif PFMs.

    Parse PFM motif files. Read data are used to build the motif PSSM,
    P-value matrix, etc. 

    ...

    Parameters
    ----------
    motif_file : str
        path motif file (PFM format)
    bg_file
        path to background file in Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html)
    pseudocount : float
        value to add to motif PWM counts
    no_reverse : bool
        if False, consider only forward strand. Use both strands otherwise
    verbose : bool
        print additional information
    debug:
        trace the full error stack

    Returns
    -------
    Motif
        Motif object 
    """
    
    motifs_raw: List = list()
    motifs: List[Motif] = list()
    counts: List[float] = list()
    alphabet: List[str] = list()
    motif_id: str = ""
    motif_name: str = ""
    motif_number: int = 0
    motif_read: int = 0

    try:
        ifstream = open(motif_file, mode="r")
        for line in ifstream:
            line = line.strip()
            if line:
                columns = line.split()
                if not columns: break  # EOF or empty file?
                fields = len(columns)
                if line.startswith("#"): continue  # comment
                elif line.startswith(">"):  # hocomoco PCM
                    if motif_number != 0 and motif_number != motif_read:
                        motifs_raw = __append_PFM(
                            motif_file,
                            motifs_raw,
                            motif_name,
                            motif_id,
                            counts,
                            alphabet,
                            debug
                        )
                        motif_read = len(motifs_raw)
                    motif_id = line[1:].strip().split()[0]
                    motif_name = motif_id 
                elif columns[0] == "Gene":
                    if motif_number != 0 and motif_number != motif_read:
                        motifs_raw = __append_PFM(
                            motif_file,
                            motifs_raw,
                            motif_name,
                            motif_id,
                            counts,
                            alphabet,
                            debug
                        )
                        motif_read = len(motifs_raw)
                    motif_name = columns[1]
                elif columns[0] == "Motif":
                    if motif_number != 0 and motif_number != motif_read:
                        motifs_raw = __append_PFM(
                            motif_file,
                            motifs_raw,
                            motif_name,
                            motif_id,
                            counts,
                            alphabet,
                            debug
                        )
                        motif_read = len(motifs_raw)
                    motif_id = columns[1]
                elif columns[0] == "Pos":
                    if fields == 5:
                        if motif_number != 0 and motif_number != motif_read:
                            motifs_raw = __append_PFM(
                                motif_file,
                                motifs_raw,
                                motif_name,
                                motif_id,
                                counts,
                                alphabet,
                                debug
                            )
                            motif_read = len(motifs_raw)
                        alphabet = columns[1:]
                        assert len(alphabet) == 4
                        if any([n not in DNA_ALPHABET for n in alphabet]):
                            errmsg = "Motif alphabet contains wrong characters.\n"
                            exception_handler(ValueError, errmsg, debug)
                elif columns[0] in DNA_ALPHABET:
                    if fields == 4:
                        alphabet = columns
                        assert len(alphabet) == 4
                        if any([n not in DNA_ALPHABET for n in alphabet]):
                            errmsg = "Motif alphabet contains wrong characters.\n"
                            exception_handler(ValueError, errmsg, debug)
                    else:
                        continue
                else:  # motif counts
                    if fields == 4:
                        count = list(map(float, columns))
                    elif fields == 5:
                        count = list(map(float, columns[1:]))
                    else: 
                        continue
                    if motif_number == motif_read:  # new motif starts here
                        motif_number += 1
                    counts.append(count)
            else:
                if motif_number != 0 and motif_number != motif_read:
                        motifs_raw = __append_PFM(
                            motif_file,
                            motifs_raw,
                            motif_name,
                            motif_id,
                            counts,
                            alphabet,
                            debug
                        )
                        motif_read = len(motifs_raw)
        if motif_number != 0 and motif_number != motif_read:
            motifs_raw = __append_PFM(
                motif_file,
                motifs_raw,
                motif_name,
                motif_id,
                counts,
                alphabet,
                debug
            )
            motif_read = len(motifs_raw)
    except:
        errmsg = "An error occurred while reading {}.\n"
        exception_handler(MotifFileReadError, errmsg.format(motif_file), debug)
    finally:
        ifstream.close()
    if os.path.isfile(bg_file): 
        bgs = readBGfile(bg_file, debug)
        bgs = pseudo_bg(bgs, no_reverse, debug)
    elif bg_file != UNIF:
        errmsg = "Unable to parse {}.\n"
        exception_handler(BGFileError, errmsg.format(bg_file), debug)
    for motif_raw in motifs_raw:
        motif_counts = pd.DataFrame(
            data=motif_raw["counts"], columns=motif_raw["alphabet"]
        ).T
        motif_width = motif_counts.shape[1]
        if bg_file == UNIF:
            bgs = get_uniformBG(motif_raw["alphabet"], debug)
            bgs = pseudo_bg(bgs, no_reverse, debug)
        motif_probs = (motif_counts / motif_counts.sum(0))
        motif_probs = norm_motif(
            motif_probs, 
            motif_width, 
            motif_raw["alphabet"], 
            debug
        )
        motif_probs = apply_pseudocount_jaspar_pfm_transfac(
            motif_counts.to_numpy(),
            motif_probs.to_numpy(),
            pseudocount,
            bgs,
            motif_width,
            motif_raw["alphabet"],
            motif_raw["nucsmap"],
            debug
        )
        motif: Motif = Motif(
            motif_probs, 
            motif_width,
            motif_raw["alphabet"],
            motif_raw["motif_id"], 
            motif_raw["motif_name"],
            motif_raw["nucsmap"]
        )
        motif.setBg(bgs) 
        motifs.append(motif)

    return motifs
        
# end of read_PFM_motif()


def __append_PFM(
    motif_file: str,
    motifs: List[dict],
    motif_name: str,
    motif_id: str,
    counts: List[float],
    alphabet: List[str], 
    debug: bool
) -> List[dict]:
    
    if not motif_name and not motif_id:
        dummy_name = os.path.abspath(motif_file).split("/")[-1]
        motif_id = dummy_name
        motif_name = dummy_name
    elif not motif_name:
        motif_name = motif_id
    elif not motif_id:
        motif_id = motif_name
    if not counts:
        errmsg = "Wrong motif read.\n"
        exception_handler(ValueError, errmsg, debug)
    if any([len(c) != len(counts[0]) for c in counts]):
        errmsg = "Motif counts width mismatch.\n"
        exception_handler(ValueError, errmsg, debug)
    if not alphabet:
        alphabet = DNA_ALPHABET
    nucsmap = dict()
    for i in range(len(alphabet)): nucsmap[alphabet[i]] = i
    motif = {
        "motif_name":motif_name,
        "motif_id":motif_id,
        "counts":counts,
        "alphabet":alphabet,
        "nucsmap":nucsmap
    }
    motifs.append(motif)
    print("Read motif %s." % (motif_name))
    return motifs

# end of __append_PFM()    


def build_TRANSFAC_motif(
    motif_file: str,
    bg_file: str,
    pseudocount: float,
    no_reverse: bool,
    verbose: bool,
    debug: bool
) -> Motif:
    """Build Motif object from TRANSFAC motif files.

    Computes PSSM scoring matrix from counts values stored in TRANSFAC 
    motif files. While computing the PSSM, the function build the 
    corresponding P-value matrix. The PSSM assigns log-likelihood scores
    to potential motif occurrences, while the P-value matrix computes 
    scores statistical significance. 

    ...

    Parameters
    ----------
    motif_file : str
        path to motif file (TRANSFAC format)
    bg_file
        path to background file in Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html).
    pseudocount : float
        value to add to motif PWM counts
    no_reverse : bool
        if False, consider only forward strand. Use both strands otherwise.
    verbose : bool
        print additional information
    debug : bool
        trace the full error stack

    Returns
    -------
    Motif
        Motif object
    """

    if not isinstance(motif_file, str):
        errmsg = "Expcted str, got {}.\n"
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
    if pseudocount <= 0:
        errmsg = "Pseudocount must be > 0.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(no_reverse, bool):
        errmsg = "Expected bool, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(no_reverse).__name__), debug)
    
    # parse motif
    motif: Motif = read_TRANSFAC_motif(
        motif_file, bg_file, pseudocount, no_reverse, verbose, debug
    )
    if verbose: start_mp = time.time()
    motif = process_motif_for_logodds(motif, debug)
    if verbose:
        stop_mp = time.time()
        msg: str = "Motif %s processed in %.2fs" % (motif.motifID, (stop_mp - start_mp))
        print(msg)
    
    return motif

# end of build build_TRANSFAC_motif()


def read_TRANSFAC_motif(
    motif_file:str,
    bg_file:str,
    pseudocount:str,
    no_reverse:bool,
    verbose:bool,
    debug:bool
) -> Motif:
    """Read motif TRANSFAC PWMs.

    Parse TRANSFAC motif files. Read data are used to build the motif PSSM,
    P-value matrix, etc. 

    ...

    Parameters
    ----------
    motif_file : str
        path motif file (TRANSFAC format)
    bg_file
        path to background file in Markov Background Format
        (http://meme-suite.org/doc/bfile-format.html)
    pseudocount : float
        value to add to motif PWM counts
    no_reverse : bool
        if False, consider only forward strand. Use both strands otherwise
    verbose : bool
        print additional information
    debug:
        trace the full error stack

    Returns
    -------
    Motif
        Motif object 
    """
  
    counts: Dict = dict()
    if verbose: start_rm = time.time()
    try:
        ifstream = open(motif_file, mode="r")
        for line in ifstream:
            line = line.strip()
            if not line: continue
            key_value = line.split(None, 1)
            key = key_value[0].strip()
            if key == "AC": 
                assert len(key_value) == 2
                motifID = key_value[1].strip()
            elif key == "ID":
                assert len(key_value) == 2
                motifName = key_value[1].strip()
            elif key in ("PO", "P0"):  # old TRANSFAC files use PO instead of P0
                assert len(key_value) == 2
                value = key_value[1].strip()
                nucs = value.split()[:4]
                assert nucs == DNA_ALPHABET
                for nuc in nucs: counts[nuc] = []
                width = 0
                for line in ifstream:
                    line = line.strip()
                    key_value = line.split(None, 1)
                    key = key_value[0].strip()
                    try:
                        pos = int(key)
                    except ValueError:
                       break  # matrix completed
                    if len(key_value) != 2:
                        errmsg = "Invalid count line in {}.\n"
                        exception_handler(
                            ValueError, errmsg.format(motif_file), debug
                        )
                    value = key_value[1].strip()
                    width += 1
                    if pos != width:
                        errmsg = "Mismatching motif position and width.\n"
                        exception_handler(ValueError, errmsg, debug)
                    values = value.split()[:4]
                    if len(values) != 4:
                        errmsg = "TRANSFAC matrix should have count values for each nucloetide.\n"
                        exception_handler(ValueError, errmsg, debug)
                    for nuc, v in zip(nucs, values): 
                        counts[nuc].append(float(v))
    except:
        errmsg = "An error occurred while reading {}.\n"
        exception_handler(MotifFileReadError, errmsg.format(motif_file), debug)
    else:
        if any([len(counts[DNA_ALPHABET[0]]) != len(counts[nuc]) for nuc in counts.keys()]):
            errmsg = "Motif counts width mismatch.\n"
            exception_handler(ValueError, errmsg, debug)
        nucsmap = dict()
        for i in range(len(nucs)): nucsmap.update({nucs[i]:i})
        motif_counts: pd.DataFrame = pd.DataFrame(counts).T
        motif_width: int = len(counts[DNA_ALPHABET[0]])
        alphabet:list = sorted(nucs)
        
        # compute background
        if bg_file == UNIF: bgs = get_uniformBG(alphabet, debug)
        elif os.path.isfile(bg_file): bgs = readBGfile(bg_file, debug)
        else:
            errmsg = "Unable to parse {}.\n"
            exception_handler(BGFileError, errmsg.format(bg_file), debug)
        bgs = pseudo_bg(bgs, no_reverse, debug)  # add pseudocount to background

        # motif probability matrix
        motif_probs = (motif_counts / motif_counts.sum(0))
        motif_probs = norm_motif(motif_probs, motif_width, alphabet, debug)
        motif_probs = apply_pseudocount_jaspar_pfm_transfac(
            motif_counts.to_numpy(), 
            motif_probs.to_numpy(),
            pseudocount,
            bgs,
            motif_width,
            alphabet,
            nucsmap,
            debug
        )                    
        motif: Motif = Motif(
            motif_probs, motif_width, alphabet, motifID, motifName, nucsmap
        )   
        motif.setBg(bgs)
        if verbose:
            stop_rm = time.time()
            msg:str = "Read motif %s in %.2fs" % (motifID, (stop_rm - start_rm))
            print(msg)
    finally:
        ifstream.close()

    return motif

# end of read_TRANSFAC_motif()


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
    motifFF: str = str()
    if not isinstance(motif_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_file).__name__), debug)
    if not os.path.isfile(motif_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(motif_file), debug)
    if isJaspar_ff(motif_file, debug):
        motifFF = JASPAR
    elif isMEME_ff(motif_file, debug):
        motifFF = MEME
    elif isPFM_ff(motif_file, debug):
        motifFF = PFM
    elif isTRANSFAC_ff(motif_file, debug):
        motifFF = TRANSFAC
    else:
        errmsg = "Motif PWM must be in MEME, JASPAR, PFM, or TRANSFAC format.\n"
        exception_handler(MotifFileFormatError, errmsg, debug)

    assert bool(motifFF)  # should have been assigned FF
    # chhose motif PWM parsing method
    if motifFF == JASPAR:
        motif = build_motif_JASPAR(
            motif_file, bgs, pseudo, no_reverse, verbose, debug
        )
    elif motifFF == MEME:
        motif = build_motif_MEME(
            motif_file, bgs, pseudo, no_reverse, cores, verbose, debug
        )
    elif motifFF == PFM:
        motif = build_PFM_motif(
            motif_file, bgs, pseudo, no_reverse, cores, verbose, debug
        )
    elif motifFF == TRANSFAC:
        motif = build_TRANSFAC_motif(
            motif_file, bgs, pseudo, no_reverse, verbose, debug
        )
    else:  # should never reach this point
        errmsg = "Unrecognized motif file format.\n"
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
        errmsg = "Unknown motif alphabet.\n"
        exception_handler(ValueError, errmsg, debug)

    # tolerance in the difference between the position probability and 1
    tolerance: float = 0.00001 
    for j in range(motif_width):
        tot = np.double(0)
        for nuc in alphabet: tot += motif_probs.loc[nuc,j]
        assert tot != 0
        if not almost_equal(1, tot, tolerance):
            for nuc in alphabet:
                motif_probs.loc[nuc,j] = np.double(motif_probs.loc[nuc,j] / tot)

    return motif_probs

# end norm_motif()

