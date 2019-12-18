"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

The script defines the MotifSet class, that stores a user defined number
of motifs which will be used to scan the genome graphs paths

"""

from . import motif
import sys
from grafimo.utils import die

class MotifSet(object):

    """
        Class to represent a set of motifs
    """

    motifs = [] # list of motifs

    def __init__(self):
        pass

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
            sys.stderr.write("ERROR: The motifs to add to the set must be in a list")
            die(1)

        if not all(isinstance(m, motif.Motif) for m in mtfs):
            sys.stderr.write("ERROR: Only an instance of Motif can be added to the motif list. Exiting...")
            die(1)

        motifs = self.motifs
        motifs += mtfs # add the motif
        self.motifs = motifs

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

        if not self.motifs:
            sys.stderr.write("ERROR: trying to get an empty list of motifs. Exiting...")
            die(1)

        return self.motifs

