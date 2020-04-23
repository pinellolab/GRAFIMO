"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

The script defines the MotifSet class, that stores a user defined number
of motifs which will be used to scan the genome graphs paths

"""

from . import motif


class MotifSet(object):
    """
        Class to represent a set of motifs
    """

    _motifs = []  # list of motifs
    _motifs_num = 0  # number of motifs in the set

    def __init__(self):
        pass

    def __iter__(self):
        return MotifSetIterator(self)

    def addMotif(self, mtfs):
        """
            Add a motif to the motif list
            ----
            Parameters:
                mtfs (list) : list of motifs to add to the motif set
            ----
            Returns:
                None
        """

        if not isinstance(mtfs, list):
            errmsg = "\n\nERROR: The motifs to add to the set must be in a list"
            raise ValueError(errmsg)

        if not all(isinstance(m, motif.Motif) for m in mtfs):
            errmsg = "\n\nERROR: Only an instance of Motif can be added to the motif list. Exiting..."
            raise ValueError(errmsg)

        motifs = self._motifs
        motifs += mtfs  # add the motif
        self._motifs = motifs

        self._motifs_num += len(mtfs)
        assert self._motifs_num > 0
    # end of addMotif()

    def getMotifsList(self):
        """
            Returns the list of motifs
            ----
            Parameters:
                None
            ----
            Returns:
                motifs (list) : list of the motifs
        """

        if not self._motifs:
            errmsg = "\n\nERROR: trying to get an empty list of motifs. Exiting..."
            raise Exception(errmsg)

        return self._motifs
    # end of getMotifsList()

    def length(self):
        assert self._motifs_num >= 0
        return self._motifs_num
# end of MotifSet


class MotifSetIterator:

    def __init__(self, motifSet):
        if not isinstance(motifSet, MotifSet):
            raise ValueError("\n\nERROR: the iterable object must be an instance of MotifSet")
        self._motifSet = motifSet  # reference to MotifSet object
        self._index = 0
    # end of __init__()

    def __next__(self):
        if self._index < (self._motifSet.length()):
            motifs_lst = self._motifSet.getMotifsList()
            result = motifs_lst[self._index]
            self._index += 1
            return result

        # end of Iteration
        raise StopIteration
    # end of __next__()
# end of MotifSetIterator

