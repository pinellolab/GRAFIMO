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
    #-------------------------------------------------------------------
    _motifs: List[Motif]
    _motifs: int

    _motifs = list()
    _motifs_num = 0 

    #-------------------------------------------------------------------
    # MotifSet methods
    #-------------------------------------------------------------------
    def __init__(self):
        pass


    def __iter__(self):
        return MotifSetIterator(self)


    def addMotif(self, mtfs: List[Motif]) -> None:
        """Add a list object containing Motif instanced to the Motifset
        object

        Parameters
        ----------
        mtfs : list
            list containing the Motif objects to add to the MotifSet
        """

        errmsg: str
        if not isinstance(mtfs, list):
            errmsg = "\n\nERROR: The motifs to add to the set must be in a list"
            raise ValueError(errmsg)

        if not all(isinstance(m, Motif) for m in mtfs):
            errmsg = "\n\nERROR: Only Motif instances can be added to MotifSet"
            raise ValueError(errmsg)

        motifs: List[Motif]
        motifs = self._motifs
        motifs += mtfs  # add the motifs
        self._motifs = motifs

        self._motifs_num += len(mtfs)  # update the motifs number
        assert self._motifs_num > 0


    def getMotifsList(self) -> List[Motif]:
        """Returns a list containing the Motif instances contained in 
        MotifSet

        Returns
        -------
        list
            list containing the Motif objects contained in MotifSet
        """

        errmsg: str
        if not self._motifs:
            errmsg = "\n\nERROR: trying to access an empty list of motifs"
            raise Exception(errmsg)

        return self._motifs
   

    def length(self) -> int:
        """Returns the number of motifs currently contained in the 
        MotifSet

        Returns
        -------
        int
            number of motifs contained innthe MotifSet
        """

        assert self._motifs_num >= 0
        return self._motifs_num

# end of MotifSet


class MotifSetIterator:
    """
    This class defines an iterator for the MotifSet class. With this 
    class a MotifSet object can be iterated with a simple for loop, 
    through the corresponding __iter()__ method.

    ...

    Attributes
    ----------
    _motifSet : MotifSet
        reference to the itereted MotifSet
    _index : int
        element index
    """

    #-------------------------------------------------------------------
    # MotifSetIterator attributes
    #-------------------------------------------------------------------
    _motifSet: MotifSet
    _index: int

    #-------------------------------------------------------------------
    # MotifSetIterator methods
    #-------------------------------------------------------------------
    def __init__(self, motifSet: MotifSet):
        errmsg: str
        if not isinstance(motifSet, MotifSet):
            errmsg = "\n\nERROR: needed an instance of MotifSet"
            raise ValueError(errmsg)

        self._motifSet = motifSet
        self._index = 0


    def __next__(self) -> Motif:
        if self._index < (self._motifSet.length()):
            motifs_lst: List[Motif] = self._motifSet.getMotifsList()
            result: Motif = motifs_lst[self._index]
            self._index += 1
            return result

        raise StopIteration
   
# end of MotifSetIterator

