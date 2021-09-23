"""Functions and constant variables used in GRAFIMO code"""

from typing import (
    List, 
    Optional, 
    Tuple, 
    Dict
)
from colorama import (
    Fore, 
    init
) 
from shutil import which

import pandas as pd
import numpy as np

import gzip
import sys
import os


#-----------------------------------------------------------------------
# constant vars
#
DNA_ALPHABET = ["A", "C", "G", "T"] 
REV_COMPL = {"A":"T", "C":"G", "G":"C", "T":"A"}
NOMAP = "NOMAP"
ALL_CHROMS = "use_all_chroms"
UNIF = "unfrm_dst"
PSEUDObg = np.double(0.0000005)
LOG_FACTOR = 1.44269504
RANGE = 1000
CHROMS_LIST = [str(i) for i in range(1, 23)] + ['X', 'Y']
DEFAULT_OUTDIR = "default_out_dir_name"
EXT_DEPS = ['tabix', 'vg', 'dot']
SOURCE = 'grafimo'
TP = 'nucleotide_motif'
PHASE = '.'
JASPAR = bin(0)
MEME = bin(1)
PFM = bin(3)
TRANSFAC = bin(4)


#-----------------------------------------------------------------------
# functions
#
def die(code: int) -> None:
    """Stop the execution and exit.

    ...

    Parameters
    ----------
    code : int
        stop code
    """

    sys.exit(code)

# end of die()


def sigint_handler() -> None:
    """Print a message when a SIGINT is caught and exit."""

    print("\nCaught SIGINT. GRAFIMO will exit")
    die(2)

# end of sigint_handler()


def exception_handler(
    exception_type, 
    exception, 
    debug
):
    """Handle Exceptions for debugging purposes. If debug is set to True 
    the entire error traceback is printed to stderr, otherwise a 
    graceful error message is sent to stderr.
    """
    init()
    if debug:
        raise exception_type("\n\n{}".format(exception))
    else:
        sys.stderr.write(Fore.RED + "\n\nERROR: " + "{}".format(exception) + Fore.RESET)
        die(1)

# end of exception_handler()


def parse_namemap(namemap_fn:str) -> Dict:
    """Parse chromosome name-map file.

    Parameters
    ----------
    namemap_fn : str
        path to name-map file

    Returns
    -------
        dict
    """

    if not isinstance(namemap_fn, str):
        errmsg = "\n\nERROR: expected str, got {}.\n"
        raise TypeError(errmsg.format(type(namemap_fn).__name__))
    chroms_namemap = dict()
    if namemap_fn == NOMAP:  # no name-map file given -> return empty dict
        assert not bool(chroms_namemap)
    else:  # name-map file is given
        if not os.path.isfile(namemap_fn):
            errmsg = "\n\nERROR: Unable to find {}.\n"
            raise FileNotFoundError(errmsg.format(namemap_fn))
        try:
            with open(namemap_fn, mode="r") as infile:
                for line in infile:
                    chrom, name = line.strip().split()
                    chroms_namemap.update({chrom:name})
        except:
            errmsg = "\n\nERROR: a problem was encountered while reading {}.\n"
            raise IOError(errmsg.format(namemap_fn))
        finally:
            infile.close()
        assert bool(chroms_namemap)
    return chroms_namemap

# end of parse_namemap()


def isListEqual(lst1: List, lst2: List) -> bool:
    """Check if two lists contain the same distinct elements.

    Parameters
    ----------
    lst1 : list
        list 1
    lst2 : list
        list 2
    
    Returns
    -------
    bool
    """

    if (len(lst1) == len(lst2) and set(lst1) == set(lst2)):
        return True

    return False

# end of isListEqual()

def anydup(lst: List) -> bool:
    """Check if are there any duplicate value in the input list
    
    Parameters
    ----------
    lst : list
        list of values

    Returns 
    -------
    bool
    """

    seen = set()
    for e in lst:
        if e in seen: return True
        seen.add(e)
    return False

# end of anydup()


def initialize_chroms_list(args_chroms: List[str]):
    """Initialize a list of chromomsome to use during GRAFIMO analysis.

    Parameters
    ----------
    args_chroms : list
        list of user chromosomes
    
    Returns
    -------
    list
        chromosome list
    """

    if not args_chroms:
        chroms = CHROMS_LIST # all the chromosomes

    else:
        chroms = args_chroms # the given chromosomes

    return chroms

# end of initialize_chroms_list()


def check_deps() -> Tuple[bool, List[str]]:
    """Check if the external dependencies needed to run GRAFIMO are 
    satisfied.

    Returns
    -------
    bool
        check result
    list
        list of dependencies not satisfied
    """

    deps_not_sats: List[str] = list()
    sat: bool = True

    for dep in EXT_DEPS:
        if not which(dep) is not None:
            deps_not_sats.append(dep)
            sat = False

    return sat, deps_not_sats

# end of check_deps()


def isJaspar_ff(motif_file: str, debug: bool) -> bool:
    """Check if the given motif file is a JASPAR file.
        
    Parameters
    ----------
    motif_file : str 
        motif file
    debug : bool
        trace the full error stack
        
    Returns
    -------
    bool
        check result
    """

    if not isinstance(motif_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_file).__name__), debug)
    if not os.path.isfile(motif_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(motif_file), debug)
    if os.stat(motif_file).st_size == 0:
        errmsg = "{} seems to be empty.\n"
        exception_handler(EOFError, errmsg.format(motif_file), debug)

    ff = motif_file.split(".")[-1]
    if ff != "jaspar":
        return False
    ifstream = open(motif_file, mode="r")
    header = ifstream.readline().strip()
    if not header.startswith(">"):  # probably not JASPAR
        return False
    while True:
        line = ifstream.readline()
        line = line.strip().split()
        if not line: break  # EOF or empty file?
        if not line[1] == "[" and line[-1] == "]":
            return False  # surely not JASPAR
        counts = line[2:-1]
        if any([not is_number(c, debug) for c in counts]):
            return False
    return True

# end of isJaspar_ff()


def isMEME_ff(motif_file: str, debug: bool) -> bool:
    """Check if the given motif file is a MEME file.
        
    Parameters
    ----------
    motif_file : str 
        motif file
    debug : bool
        trace the full error stack
        
    Returns
    -------
    bool
        check result
    """

    if not isinstance(motif_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_file).__name__), debug)
    if not os.path.isfile(motif_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(motif_file), debug)
    if os.stat(motif_file).st_size == 0:
        errmsg = "{} seems to be empty.\n"
        exception_handler(EOFError, errmsg.format(motif_file), debug)

    ifstream = open(motif_file, mode="r")
    for line in ifstream:
        if line.startswith("MEME version"): return True
    return False  # no MEME version found --> improper input

# end of isMEME_ff()


def isPFM_ff(motif_file: str, debug: bool) -> bool:
    """Check wheter motif files are in PFM format.
    Assumes that each file contains at most one motif.

    ...

    Parameters
    ----------
    motif_file: str
        motif file
    debug: bool
        print the full error stack
    
    Returns
    -------
    bool
        check result
    """

    if not isinstance(motif_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_file).__name), debug)
    if not os.path.isfile(motif_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(motif_file), debug)
    if os.stat(motif_file).st_size == 0:
        errmsg = "{} seems to be empty.\n"
        exception_handler(EOFError, errmsg.format(motif_file), debug)

    ifstream = open(motif_file, mode="r")
    for line in ifstream:
        if line.startswith(">"): continue  # probably from JASPAR
        counts = line.strip().split()
        if any([not is_number(c, debug) for c in counts]):
            return False
    return True


def isTRANSFAC_ff(motif_file: str, debug: bool) -> bool:
    """Check wheter motif files are in TRANSFAC format.
    Assumes that each file contains at most one motif.

    ...

    Parameters
    ----------
    motif_file: str
        motif file
    debug:
        print the full error stack

    Returns
    -------
    bool
        check result
    """

    if not isinstance(motif_file, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(motif_file).__name__), debug)
    if not os.path.isfile(motif_file):
        errmsg = "Unable to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg.format(motif_file), debug)
    if os.stat(motif_file).st_size == 0:
        errmsg = "{} seems to be empty.\n"
        exception_handler(EOFError, errmsg.format(motif_file), debug)

    transfac_key_fields = {
        "AC":False,
        "ID":False,
        "PO":False,

    }
    width = 0
    ifstream = open(motif_file, mode="r")
    for line in ifstream:
        line = line.strip()
        if not line: continue  # empty lines allowed
        line = line.split(None, 1)
        key = line[0].strip()
        if len(key) != 2:
            return False
        if len(line) == 2:  # key - value
            value = line[1].strip()
            if key in transfac_key_fields.keys():
                if not value: return False
                if key in ("P0", "PO"):  # old TRANSFAC files use PO instead of P0
                    motif_alphabet = value.split()[:4]
                    if motif_alphabet != DNA_ALPHABET: return False
                transfac_key_fields[key] = True
            try:
                pos = int(key)  # check if counts line
            except:
                continue  
            else:
                if width == 0 and pos == 0:
                    return False  # TRANSFAC motif should start from 0
                width += 1
                if width != pos: return False  # position mismatch
    if sum(transfac_key_fields.values()) == 3: return True
    return False
            
            



    


def isbed(bedfile: str, debug: bool) -> bool:
    """Check if the given file is in UCSC BED format.

    ...

    Parameters
    ----------
    bedfile : str
        path to bedfile
    debug : bool
        trace the full error stack

    Returns
    -------
    bool
        check result
    """

    if not isinstance(bedfile, str):
        errmsg = "Expected str, got {}.\n"
        exception_handler(TypeError, errmsg, debug)
    if not os.path.isfile(bedfile):
        errmsg = "Unble to locate {}.\n"
        exception_handler(FileNotFoundError, errmsg, debug)
    
    # check if bed is gzipped
    gzipped: bool = False
    ff: str = bedfile.split(".")[-1]
    if ff == "gz": gzipped = True
    if gzipped: ifstream = gzip.open(bedfile, mode="rt")
    else: ifstream = open(bedfile, mode="r")
    for line in ifstream:
        if line.startswith("chr"): 
            if len(line.split()) >= 3: return True  # at least chrom, start, end
            else: return False
    # if EOF reached and no line started with chr --> False
    return False


def almost_equal(
    value1: np.double,
    value2: np.double,
    slope: np.double
) -> bool:
    """Check if two values are close to each other, given a degree
    of tolerance (slope).

    ...
    
    Parameters
    ----------
    value1 : np.double
        first value
    value2 : np.double 
        second value
    slope : np.double
        tolerance
        
    Returns
    -------
    bool
    """

    if ((value1 - slope) > value2 or (value1 + slope) < value2):
        return False
    return True

# end of almost_equal()


def lg2(value: np.double) -> np.double:
    """C like computation of log2
        
    Parameters
    ----------
    value : np.double 
        value
        
    Returns
    np.double
        log2 of the input value
    """

    log2value = (np.log(value) * LOG_FACTOR)
    return log2value

# end of lg2()


def dftolist(data: pd.DataFrame, no_qvalue: bool, debug: bool) -> List:
    """Convert pandas DataFrame in a list of lists.

    ...

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame
    no_qvalue : bool
        skip q-values
    debug : bool
        trace the full error stack

    Returns
    -------
    list
        DataFrame values as list of lists
    """

    if not isinstance(data, pd.DataFrame):
        errmsg = "Expected pandas.DataFrame, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(data).__name__), debug)
    if len(data) == 0:
        errmsg = "Empty DataFrames cannot be converted to lists of values.\n"
        exception_handler(ValueError, errmsg, debug)
    if len(data.columns) > 12 or len(data.columns) < 11:
        errmsg = "Not enough values to extract from the DataFrame.\n"
        exception_handler(ValueError, errmsg, debug)
    if not isinstance(no_qvalue, bool):
        errmsg = "Expected bool, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(no_qvalue).__name__), debug)

    seqnames: List[str] = data["sequence_name"].tolist()
    starts: List[int] = data["start"].tolist()
    stops: List[int] = data["stop"].tolist()
    scores: List[np.double] = data["score"].tolist()
    strands: List[str] = data["strand"].tolist()
    motifIDs: List[str] = data["motif_id"].tolist()
    motifNames: List[str] = data["motif_alt_id"].tolist()
    pvalues: List[np.double] = data["p-value"].tolist()
    sequences: List[str] = data["matched_sequence"].tolist()
    frequencies:List[int] = data["haplotype_frequency"].tolist()
    references: List[str] = data["reference"].tolist()
    if not no_qvalue:
        qvalues: List[np.double] = data["q-value"].tolist()
        summary = [
            motifIDs, motifNames, seqnames, starts, stops, strands, scores,
            pvalues, sequences, frequencies, references, qvalues
        ]
    else:  # no_qvalue == True
        summary = [
            motifIDs, motifNames, seqnames, starts, stops, strands, scores,
            pvalues, sequences, frequencies, references
        ]
    if any([len(seqnames) != len(l) for l in summary]):
        errmsg = "List length mismatch.\n"
        exception_handler(ValueError, errmsg, debug)

    return summary

# end of dftolist()


def is_number(s: str, debug: bool) -> bool:
    """Check wheter a string contains float values.

    ...

    Parameters
    ----------
    s : str
        string object
    debug : bool
        trace the full error stack
    
    Returns
    -------
    bool
    """

    if not isinstance(s, str):
        errmsg = "\n\nExpected str, got {}.\n"
        exception_handler(TypeError, errmsg.format(type(s).__name__), debug)
    
    try: 
        float(s)
        return True
    except ValueError:
        return False

# end of is_number()


def printProgressBar(
    iteration: int,
    total: int,
    prefix: Optional[str] = '',
    suffix: Optional[str] = '',
    decimals: Optional[int] = 1,
    length: Optional[int] = 50,
    fill: Optional[str] = '=',
    printEnd: Optional[str] ="\r"
) -> None:
    """Print progress bar.

    ...

    Parameters
    ----------
    iteration : int
        iteration number
    total : int
        total number of iterations 
    prefix : str
        progress bar prefix
    suffix : str
        progress bar suffix
    decimals : int
        number of decimal digits to display
    length : int
        bar length
    fill : str
        bar filling character
    printEnd : str
        bar ending character 
    """

    pct: str = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength: int = int(length * iteration // total)
    # allocate space for the bar
    bar: float = fill * filledLength + ' ' * (length - filledLength)  
    print('\r%s [%s] %s%% %s' % (prefix, bar, pct, suffix), end = printEnd)
    # new line when the bar is completely filled
    if iteration == total:
        print()

# end of printProgressBar()

