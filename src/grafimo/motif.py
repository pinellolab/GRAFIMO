"""Motif object definition.

The motif PWM in JASPAR or MEME format is represented by the motif class
in GRAFIMO. In a single motif object can be accessed the corresponding
probability matrices, scaled score matrices, P-value matrices, etc.
"""


from grafimo.GRAFIMOException import NotValidMotifMatrixException, \
    NoDataFrameException, WrongMotifIDException, WrongMotifWidthException, \
    WrongMotifNameException, NotValidAlphabetException, ValueException, \
    NotValidBGException
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
    _count_matrix : pandas.DataFrame
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

    #-------------------------------------------------------------------
    # Motif attributes
    # ------------------------------------------------------------------
    _count_matrix: pd.DataFrame  
    _score_matrix: np.ndarray
    _pval_matrix: np.array
    _min_val: np.double
    _max_val: np.double
    _scale: int
    _offset: np.double
    _bg: dict  
    _width: int
    _motif_id: str 
    _motif_name: str  
    _alphabet: list
    _isScaled: bool

    # class attributes value initialization
    _min_val = -np.inf
    _max_val = np.inf
    _scale = -1
    _offset = 0
    _width = -1
    _isScaled = False

    #-------------------------------------------------------------------
    # Motif methods
    # ------------------------------------------------------------------
    def __init__(self,
                 count_matrix: pd.DataFrame,
                 width: int,
                 alphabet: List[str],
                 motif_id: str,
                 motif_name: str
    ):

        errmsg: str

        if count_matrix.empty:
            errmsg = "\n\nERROR: attempt to initialize the motif object with an "
            errmsg += "empty count matrix"
            raise NotValidMotifMatrixException(errmsg)

        if not isinstance(count_matrix, pd.DataFrame):
            errmsg = "\n\nERROR: the given value is not a pandas.DatFrame instance"
            raise NoDataFrameException(errmsg)

        if not isinstance(width, int) or width < 0:
            errmsg = "\n\nERROR: attempt to initialize motif without a valid width"
            raise WrongMotifWidthException(errmsg)

        if not isinstance(motif_id, str) or not motif_id:
            errmsg = "\n\nERROR: cannot initialize the motif with the given ID"
            raise WrongMotifIDException(errmsg)

        if not isinstance(motif_name, str) or not motif_name:
            errmsg = "\n\nERROR: cannot initialize the motif with the given name"
            raise WrongMotifNameException(errmsg)

        if not isinstance(alphabet, list) or not isListEqual(alphabet, DNA_ALPHABET):
            errmsg = "\n\nERROR: cannot initialize a motif object with a wrong alphabet"
            raise NotValidAlphabetException(errmsg)

        self._count_matrix = count_matrix
        self._width = width
        self._motif_id = motif_id
        self._motif_name = motif_name
        self._alphabet = alphabet


    def setMotif_matrix(self,
                        motif_matrix: pd.DataFrame
    ) -> None:

        errmsg: str

        if motif_matrix.empty:
            errmsg = "\n\nERROR: attempt to use an empty motif matrix"
            raise NotValidMotifMatrixException(errmsg)

        if not isinstance(motif_matrix, pd.DataFrame):
            errmsg = "\n\nERROR: the given value is not a pandas.DataFrame instance"
            raise NoDataFrameException(errmsg)

        self._count_matrix = motif_matrix


    def setMotif_scoreMatrix(self,
                             score_matrix: np.ndarray
    ) -> None:
    
        errmsg: str
        if (not isinstance(score_matrix, np.ndarray) and
                not isinstance(score_matrix, pd.DataFrame)):
            errmsg = "\n\nERROR: the given data-structure is not an instance of "
            errmsg += "numpy.ndarray or pandas.DataFrame"
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
                             pval_mat: np.array
    ) -> None:
       
        # empty or not valid p-value matrix
        if len(pval_mat) == 0 or sum(pval_mat[:]) <= 0:
            errmsg = "\n\nERROR: invalid p-value matrix"
            raise NotValidMotifMatrixException(errmsg)

        self._pval_matrix = pval_mat


    def setMin_val(self,
                   min_val: int
    ) -> None:
        
        if min_val <= -np.inf:
            errmsg = ' '.join(["\n\nERROR: impossible to assign", min_val, 
                               "to Motif.min_val"])
            raise ValueException(errmsg)

        self._min_val = min_val


    def setMax_val(self,
                   max_val: int
    ) -> None:

        if max_val >= np.inf:
            errmsg = ' '.join(["\n\nERROR: impossible to assign", max_val, 
                               "to Motif.max_val"])
            raise ValueException(errmsg)

        self._max_val = max_val


    def setScale(self,
                 scale: int
    ) -> None:

        if not isinstance(scale, int):
            raise ValueException("\n\nERROR: the scale factor must be an int")

        assert scale > 0

        self._scale = scale


    def setOffset(self,
                  offset: np.double
    ) -> None:    
        self._offset = offset


    def setBg(self,
              bgs: Dict
    ) -> None:
        
        if not isinstance(bgs, dict):
            errmsg = "\n\nERROR: the background values are not in a dictionary"
            raise NotValidBGException(errmsg)

        self._bg = bgs


    def setWidth(self,
                 width: int
    ) -> None:
        
        if not isinstance(width, int) or width <= 0:
            errmsg = "\n\nERROR: attempt to initialize motif without a valid width"
            raise WrongMotifWidthException(errmsg)

        self._width = width


    def setMotifID(self,
                   motif_id: str
    ) -> None:
    
        if not isinstance(motif_id, str):
            errmsg = ' '.join(["\n\nERROR: cannot initialize the motif with the given ID:", 
                               motif_id])
            raise WrongMotifIDException(errmsg)

        if not motif_id:
            errmsg = "\n\nERROR: cannot use a motif with an empty ID"
            raise WrongMotifIDException(errmsg)

        self._motif_id = motif_id

    def setMotifName(self,
                     motif_name
    ) -> None:
        
        if not isinstance(motif_name, str):
            errmsg = ' '.join(["Cannot initialize the motif with the given name:", 
                               motif_name])
            raise WrongMotifNameException(errmsg)

        if not motif_name:
            errmsg = "\n\nERROR: cannot use a motif with an empty name"
            raise WrongMotifNameException(errmsg)

        self._motif_name = motif_name


    def setAlphabet(self,
                    alphabet: List[str]
    ) -> None:

        if not isinstance(alphabet, list):
            errmsg = "\n\nERROR: the given alphabet is not in a list"
            raise NotValidAlphabetException(errmsg)

        if not isListEqual(alphabet, DNA_ALPHABET):
            errmsg = "\n\nERROR: the given alphabet is not a valid DNA alphabet"
            raise NotValidAlphabetException(errmsg)

        self.alphabet = alphabet


    def setIsScaled(self,
                    isScaled: bool
    ) -> None:
        
        if not isinstance(isScaled, bool):
            raise Exception("\n\nERROR: the isScaled value must be a boolean")

        self._isScaled = isScaled


    def getMotif_matrix(self) -> pd.DataFrame:
        return self._count_matrix


    def getMotif_scoreMatrix(self) -> np.ndarray:
        return self._score_matrix


    def getMotif_pval_mat(self) -> np.array:
        return self._pval_matrix


    def getMin_val(self) -> int:
        return self._min_val


    def getMax_val(self) -> int:
        return self._max_val


    def getScale(self) -> int:
        return self._scale


    def getOffset(self) -> np.double:
        return self._offset


    def getBg(self) -> dict:
        return self._bg


    def getWidth(self) -> int:
        return self._width


    def getMotifID(self) -> str:
        return self._motif_id


    def getMotifName(self) -> str:
        return self._motif_name


    def getAlphabet(self) -> List[str]:
        return self._alphabet


    def getIsScaled(self) -> bool:
        return self._isScaled


    def compute_minValue(self) -> None:
        motif_matrix = self.getMotif_matrix()
        min_value = motif_matrix.min().sum()
        self._min_val = min_value


    def print(self,
              matrix: str
    ) -> None:
        
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
    
# end of Motif

