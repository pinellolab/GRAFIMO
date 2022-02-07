"""Motif object definition.

A single Motif object stores important values, used during the following steps of 
GRAFIMO's analysis, such as motif PSSM, p-value matrix, scaling factor, offset,
motif information, etc.
"""


from grafimo.grafimo_errors import NotValidMotifMatrixError 
from grafimo.utils import isListEqual, DNA_ALPHABET  

from typing import List, Optional, Dict   

import pandas as pd
import numpy as np


class Motif(object):
    """
    This class defines a DNA motif object.

    In a single object we carry: 
    * the original count matrix or probability matrix 
    * the motif scaled scoring matrix 
    * the P-value matrix used to assign a P-value to each motif 
      occurrence candidate score 
    * the parameters used to scale the matrix (to revert the scaled 
      score to the log-odds score) 
    * the background probability distribution used, while processing the
      PWM values 
    * the motif width
    * the minimum value in the scoring matrix
    * the maximum value in the scoring matrix
    * the motif name (both ID and extended name)
    * the motif alphabet 
        
    ...

    Attributes
    ----------
    _count_matrix : numpy.ndarray
        motif probability matrix
    _score_matrix : numpy.ndarray
        scaled motif scoring matrix
    _min_val : int
        minimum value of the scaled scoring matrix
    _max_value : int
        maximum value of the scaled scoring matrix
    _scale : int
        scaling value
    _offset : numpy.double
        offset used during motif matrix scaling
    _bg : dict
        background probability distribution
    _width : int
        motif width
    _motif_id : str
        motif ID
    _motif_name : str
        motif extended name
    _alphabet : list()
        DNA motif alphabet
    _isScaled : bool
        flag value to state if the scoring matrix has been scaled

    Methods
    -------
    setMotif_matrix(motif_matrix : pandas.DataFrame)
        set the count matrix
    setMotif_scoreMatrix(score_matrix : numpy.ndarray)
        set the scoring matrix
    setMotif_pval_matrix(pval_mat : numpy.array)
        set the P-value matrix
    setMin_val(min_val : int)
        set the scoring matrix minimum value
    setMax_val(max_val : int)
        set the scoring matrix maximum value
    setScale(scale : int)
        set the scoring matrix scaling factor
    setOffset(offset : numpy.double)
        set the scaling offset
    setBg(bgs : dict)
        set the background probability distribution
    setWidth(width : int)
        set motif width
    setMotifID(motif_id : str)
        set motif ID
    setMotifName(motif_name : str)
        set motif extended name
    setAlphabet(alphabet : list)
        set DNA motif alphabet
    setIsScaled(isScaled : bool)
        set the isScaled flag value
    getMotif_matrix()
        return the motif count matrix
    getMotif_scoreMatrix()
        return the motif scaled scoring matrix
    getMotif_pval_mat()
        return the P-value matrix
    getMin_val()
        return the scoring matrix minimum value
    getMax_val()
        return the scoring matrix maximum value
    getScale()
        return the matrix scaling factor
    getOffset()
        return the offset used while scaling the motif scoring matrix
    getBg()
        return the background probability distribution
    getWidth():
        return motif width
    getMotifID()
        return the motif ID
    getMotifName()
        return the motif extended name
    getAlphabet()
        return the DNA motif alphabet
    getIsScaled()
        return the isScaled flag value
    compute_minValue()
        compute the minimum value of the scaled scoring motif matrix
    print()
        print one matrix among the counts one, the scoring one or the 
        P-value one  
    """

    # class attributes value initialization
    _min_val = -np.inf
    _max_val = np.inf
    _scale = -1
    _offset = 0
    _width = -1
    _is_scaled = False
    #-------------------------------------------------------------------
    # Motif methods
    # 
    # these errors should never appear --> no need for error formatting
    # can assume that debug mode == True
    def __init__(
        self,
        count_matrix: np.ndarray,
        width: int,
        alphabet: List[str],
        motif_id: str,
        motif_name: str,
        nucsmap: dict
    ):
        if not isinstance(count_matrix, np.ndarray):
            errmsg = f"\n\nERROR: Expected {type(np.ndarray).__name__}, got {type(count_matrix).__name__}.\n"
            raise TypeError(errmsg)
        if count_matrix.size == 0 or sum(sum(count_matrix)) == 0:
            errmsg = "\n\nERROR: Empty motif count matrix.\n"
            raise NotValidMotifMatrixError(errmsg)
        if not isinstance(width, int):
            errmsg = f"\n\nERROR: Expected {int.__name__}, got {type(width).__name__}.\n"
            raise TypeError(errmsg)
        if width <= 0:
            errmsg = f"\n\nERROR: Forbidden motif width ({width}).\n"
            raise ValueError(errmsg)
        if not isinstance(motif_id, str):
            errmsg = f"\n\nERROR: Expected {str.__name__}, got {type(motif_id).__name__}.\n"
            raise TypeError(errmsg)
        if not motif_id:
            errmsg = "\n\nERROR: Not valid motif ID.\n"
            raise ValueError(errmsg)
        if not isinstance(motif_name, str):
            errmsg = f"\n\nERROR: Expected {str.__name__}, got {type(motif_name).__name__}.\n"
            raise TypeError(errmsg)
        if not motif_name:
            errmsg = "\n\nERROR: Not valid motif name.\n"
            raise ValueError(errmsg)
        if not isinstance(alphabet, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(alphabet).__name__}.\n"
            raise TypeError(errmsg)
        if not isListEqual(alphabet, DNA_ALPHABET):
            errmsg = "\n\nERROR: The motif is not built on DNA alphabet.\n"
            raise ValueError(errmsg)
        if not isinstance(nucsmap, dict):
            errmsg = f"\n\nERROR: Expected {dict.__name__}, got {type(nucsmap).__name__}.\n"
            raise TypeError(errmsg)
        self._count_matrix = count_matrix
        self._width = width
        self._motif_id = motif_id
        self._motif_name = motif_name
        self._alphabet = alphabet
        self._nucsmap = nucsmap


    def set_motif_matrix(self, motif_matrix: pd.DataFrame) -> None:
        if not isinstance(motif_matrix, pd.DataFrame):
            errmsg = f"\n\nERROR: Expected {type(pd.DataFrame).__name__}, got {type(motif_matrix).__name__}.\n"
            raise TypeError(errmsg)
        if motif_matrix.empty:
            errmsg = "\n\nERROR: Empty motif matrix.\n"
            raise ValueError(errmsg)
        self._count_matrix = motif_matrix


    def set_motif_score_matrix(self, score_matrix: np.ndarray) -> None:
        if not isinstance(score_matrix, np.ndarray):
            errmsg = f"\n\nERROR: Expected {type(np.ndarray).__name__}, got {type(score_matrix).__name__}.\n"
            raise TypeError(errmsg)
        if score_matrix.size == 0 or sum(sum(score_matrix)) == 0:
            errmsg = "\n\nERROR: Empty motif score matrix.\n"
            raise ValueError(errmsg)
        self._score_matrix = score_matrix


    def set_motif_pval_matrix(self, pval_mat: np.array) -> None:
        if not isinstance(pval_mat, np.ndarray):
            errmsg = f"\n\nERROR: Expected {type(np.array).__name__}, got {type(pval_mat).__name__}.\n"
            raise TypeError(errmsg)
        if len(pval_mat) == 0:
            errmsg = "\n\nERROR: Empty motif p-value matrix.\n"
            raise ValueError(errmsg)
        if sum(pval_mat) == 0:
            errmsg = "\n\nERROR: Not valid motif p-value matrix.\n"
            raise ValueError(errmsg)
        self._pval_matrix = pval_mat


    def set_min_val(self, min_val: int) -> None:
        if not isinstance(min_val, int):
            errmsg = f"\n\nERROR: Expected {int.__name__}, got {type(min_val).__name__}.\n"
            raise TypeError(errmsg)
        if min_val <= -np.inf:
            errmsg = f"\n\nERROR: Forbidden value {min_val}.\n"
            raise ValueError(errmsg)
        self._min_val = min_val


    def set_max_val(self, max_val: int) -> None:
        if not isinstance(max_val, int):
            errmsg = f"\n\nERROR: Expected {int.__name__}, got {type(max_val).__name__}.\n"
            raise TypeError(errmsg)
        if max_val >= np.inf:
            errmsg = f"\n\nERROR: Forbidden value ({max_val}).\n"
            raise ValueError(errmsg)
        self._max_val = max_val


    def set_scale(self, scale: int) -> None:
        if not isinstance(scale, int):
            errmsg = f"\n\nERROR: Expected {int.__name__}, got {type(scale).__name__}.\n"
            raise TypeError(errmsg)
        if scale <= 0:
            errmsg = "\n\nERROR: Scaling factor must be positive integer number.\n"
            raise ValueError(errmsg)
        self._scale = scale


    def set_offset(self, offset: np.double) -> None:  
        if not isinstance(offset, np.double):
            errmsg = f"\n\nERROR: Expected {type(np.double).__name__}, got {type(offset).__name__}.\n"
            raise TypeError(errmsg)  
        self._offset = offset


    def set_bg(self, bgs: Dict[str, float]) -> None:
        if not isinstance(bgs, dict):
            errmsg = f"\n\nERROR: Expected {dict.__name__}, got {type(bgs).__name__}.\n"
            raise TypeError(errmsg)
        self._bg = bgs


    def set_width(self, width: int) -> None:
        if not isinstance(width, int):
            errmsg = f"\n\nERROR: Expected {int.__name__}, got {type(width).__name__}.\n"
            raise TypeError(errmsg)
        if width <= 0:
            errmsg = "\n\nERROR: Not valid motif width.\n"
            raise ValueError(errmsg)
        self._width = width


    def set_motif_id(self, motif_id: str) -> None:
        if not isinstance(motif_id, str):
            errmsg = f"\n\nERROR: Expected {str.__name__}, got {type(motif_id).__name__}.\n"
            raise TypeError(errmsg)
        if not motif_id:
            errmsg = "\n\nERROR: Not valid motif ID.\n"
            raise ValueError(errmsg)
        self._motif_id = motif_id


    def set_motif_name(self, motif_name: str) -> None:
        if not isinstance(motif_name, str):
            errmsg = f"\n\nERROR: Expected {str.__name__}, got {type(motif_name).__name__}.\n"
            raise TypeError(errmsg)
        if not motif_name:
            errmsg = "\n\nERROR: Not valid motif name.\n"
            raise ValueError(errmsg)
        self._motif_name = motif_name


    def set_alphabet(self, alphabet: List[str]) -> None:
        if not isinstance(alphabet, list):
            errmsg = f"\n\nERROR: Expected {list.__name__}, got {type(alphabet).__name__}.\n"
            raise TypeError(errmsg)
        if len(alphabet) == 0:
            errmsg = "\n\nERROR: Empty motif alphabet.\n"
            raise ValueError(errmsg)
        if not isListEqual(alphabet, DNA_ALPHABET):
            errmsg = "\n\nERROR: The motif is not built on DNA alphabet.\n"
            raise ValueError(errmsg)
        self.alphabet = alphabet


    def set_is_scaled(self) -> None:
        if self._is_scaled:
            errmsg = "\n\nERROR: The motif matrix has already been scaled.\n"
            raise AssertionError(errmsg)
        self._is_scaled = True


    def _get_motif_matrix(self) -> np.ndarray:
        if self._count_matrix.size == 0 or sum(sum(self._count_matrix)) == 0:
            errmsg = "\n\nERROR: \"self._count_matrix\" is empty.\n"
            raise AttributeError(errmsg)
        else:
            return self._count_matrix
    
    @property
    def count_matrix(self):
        return self._get_motif_matrix()
    

    def _get_motif_score_matrix(self) -> np.ndarray:
        if self._score_matrix.size == 0 or sum(sum(self._score_matrix)) == 0:
            errmsg = "\n\nERROR: \"self._score_matrix\" is empty.\n"
            raise AttributeError(errmsg)
        else:
            return self._score_matrix

    @property
    def score_matrix(self):
        return self._get_motif_score_matrix()


    def _get_motif_pval_mat(self) -> np.ndarray:
        if self._pval_matrix.size == 0 or sum(self._pval_matrix) == 0:
            errmsg = "\n\nERROR: \"self._pval_matrix\" is empty.\n"
            raise AttributeError(errmsg)
        else:
            return self._pval_matrix
    
    @property
    def pval_matrix(self):
        return self._get_motif_pval_mat()


    def _get_min_val(self) -> int:
        return self._min_val

    @property
    def min_val(self):
        return self._get_min_val()


    def _get_max_val(self) -> int:
        return self._max_val

    @property
    def max_val(self):
        return self._get_max_val()


    def _get_scale(self) -> int:
        return self._scale

    @property
    def scale(self):
        return self._get_scale()


    def _get_nucsmap(self):
        if not bool(self._nucsmap):
            errmsg = "\n\nERROR: \"self._nucsmap\" is empty.\n"
            raise AttributeError(errmsg)
        else:
            return self._nucsmap
    
    @property
    def nucsmap(self):
        return self._get_nucsmap()


    def _get_offset(self) -> np.double:
        return self._offset

    @property
    def offset(self):
        return self._get_offset()


    def _get_bg(self) -> dict:
        if not bool(self._bg):
            errmsg = "\n\nERROR: \"self._bg\" is empty.\n"
            raise AttributeError(errmsg)
        else:
            return self._bg

    @property
    def bg(self):
        return self._get_bg()


    def _get_width(self) -> int:
        return self._width

    @property
    def width(self):
        return self._get_width()


    def _get_motif_id(self) -> str:
        if not self._motif_id:
            errmsg = "\n\nERROR: \"self._motif_id\" is empty.\n"
            raise AttributeError(errmsg)
        else:
            return self._motif_id

    @property
    def motif_id(self):
        return self._get_motif_id()


    def _get_motif_name(self) -> str:
        if not self._motif_name:
            errmsg = "\n\nERROR: \"self._motif_name\" is empty.\n"
            raise AttributeError(errmsg)
        else:
            return self._motif_name

    @property
    def motif_name(self):
        return self._get_motif_name()


    def _get_alphabet(self) -> List[str]:
        if not self._alphabet:
            errmsg = "\n\nERROR: \"self._alphabet\" si empty.\n"
            raise AttributeError(errmsg)
        else:
            return self._alphabet
    
    @property
    def alphabet(self):
        return self._get_alphabet()


    def _get_is_scaled(self) -> bool:
        return self._is_scaled

    @property
    def is_scaled(self):
        return self._get_is_scaled()


    def compute_min_value(self) -> None:
        min_value = self._score_matrix.min()
        self._min_val = min_value


    def print(self, matrix: str) -> None:
        if not isinstance(matrix, str):
            errmsg = "\n\nERROR: Expected str, got {}.\n"
            raise TypeError(errmsg.format(type(matrix).__name__))
        if not matrix:
            errmsg = "\n\nERROR: Unable to guess what should be printed.\n"
            raise ValueError(errmsg)
        available_matrices = ["raw_counts", "score_matrix", "pval_matrix"]
        if matrix not in available_matrices:
            errmsg = "\n\nERROR: Unknown motif matrix.\n"
            raise ValueError(errmsg)
        if matrix == "raw_counts": print(self._count_matrix)
        elif matrix == "score_matrix": print(self._score_matrix)
        elif matrix == "pval_matrix": print(self._pval_matrix)
        else:  # we should never reach this point
            errmsg = "\n\nERROR: Unknown motif matrix.\n"
            raise ValueError(errmsg)
    
# end of Motif

