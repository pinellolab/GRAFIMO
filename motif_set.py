"""MotifSet class and MotifSetIterator class definiton. 

A MotifSet object contains a set of DNA motifs whose occurrences will be
serached on the given genome variation graph.

This class is useful to manage the case in which the user wants to 
find the occurrences of more than one single motif.
"""

from grafimo.motif import Motif

from typing import List


class MotifSet(object):
    """
    This class represents a set of DNA motifs given as PWMs.

    A MotifSet object is an iterabole container for Motif objects. It is 
    particularly useful when searching with a single set of genomic 
    coordinates more than just a single motif.

    In the case a MEME file conatining data of more than one motif, a
    MotifSet object will contain the data for each of the available
    DNA motf, storing them as different Motif objects.

    ...

    Attributes
    ----------
    _motifs : list
        list of motifs conatined in the current object
    _motifs_num : int
        number of Motif objects contained in the current object

    Methods
    -------
    addMotif(mtfs : List[Motif])
        add the motifs in the list to the MotifSet
    getMotifsList()
        returns a list object containing all the motifs contained in
        the MotifSet
    length()
        returns the length of the MotifSet or the number of motifs 
        currently contained in the Motifset 
    """

    #-------------------------------------------------------------------
    # MotifSet attributes
    #

    _motifs = list()
    _motifs_num = 0 

    #-------------------------------------------------------------------
    # MotifSet methods
    #
    def __init__(self):
        pass


    def __iter__(self):
        return MotifSetIterator(self)

    
    def __len__(self):
        return len(self._motifs)


    def addMotif(self, mtfs: List[Motif]) -> None:
        if not isinstance(mtfs, list):
            errmsg = "\n\nERROR: Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(mtfs).__name__))
        if not all(isinstance(m, Motif) for m in mtfs):
            errmsg = "\n\nERROR: One of the list elements is not of Motif type.\n"
            raise ValueError(errmsg)
        self._motifs += mtfs
        self._motifs_num += len(self._motifs)  # update the motifs number
        assert self._motifs_num > 0


    def _get_motifs(self) -> List[Motif]:
        if not self._motifs:
            errmsg = "\n\nERROR: \"self._motifs\" is empty.\n"
            raise AttributeError(errmsg)
        return self._motifs

    @property
    def motifs(self):
        return self._get_motifs()
   

    def _get_size(self) -> int:
        assert self._motifs_num >= 0
        return self._motifs_num

    @property
    def size(self):
        return self._get_size()

# end of MotifSet


class MotifSetIterator:
    """Iterator definition for MotifSet objects. 

    ...

    Attributes
    ----------
    _motifSet : MotifSet
        MotifSet object
    _index : int
        index
    """

    #-------------------------------------------------------------------
    # MotifSetIterator attributes
    #
    _motifSet: MotifSet
    _index: int

    #-------------------------------------------------------------------
    # MotifSetIterator methods
    #
    def __init__(self, motifSet: MotifSet):
        if not isinstance(motifSet, MotifSet):
            errmsg = "\n\nERROR: Expected MotifSet, got {}.\n"
            raise TypeError(errmsg.format(type(motifSet).__name__))
        self._motifSet = motifSet
        self._index = 0


    def __next__(self) -> Motif:
        if self._index < (self._motifSet.size):
            result = self._motifSet.motifs[self._index]
            self._index += 1
            return result
        raise StopIteration
   
# end of MotifSetIterator

