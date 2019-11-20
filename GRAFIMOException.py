"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

The script contains definitions of custom Exceptions
for GRAFIMO

"""

class GRAFIMOException(Exception):
    """
        Base class for other exceptions
    """

    pass

class NoDataFrameException(GRAFIMOException):
    """
        Raise when a given object is not an instance of pandas.DataFrame
    """
    def __init__(self, message):

        self.message = message

    def __str__(self):

        return repr(self.message)

class WrongMotifWidthException(GRAFIMOException):
    """
        Raise an exception when the motif is initialized with a
        not valid width
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        return repr(self.message)

class WrongMotifIDException(GRAFIMOException):
    """
        Raise an exception when the motif object is initialized
        with a wrong ID
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        return repr(self.message)

class WrongMotifNameException(GRAFIMOException):
    """
        Raise an exception when the motif object is initialized
        with a wrong name
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        return repr(self.message)

class NotValidMotifMatrixException(GRAFIMOException):
    """
        Raise an exception when is given a not valid motif matrix
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        return repr(self.message)

class NotValidBGException(GRAFIMOException):
    """
        Raise  an exception when a not valid background distribution is
        given to the motif object
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        return repr(self.message)

class NotValidAlphabetException(GRAFIMOException):
    """
        Raise an exception when a not valid alphabet is given
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        return repr(self.message)

class NotValidFFException(GRAFIMOException):
    """
        Raise an exception if it is given a file in a not allowed
        file format
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        return repr(self.message)

class FileReadingException(Exception):
    """
        Raise an exception if something went wrong during a file
        reading
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        return repr(self.message)

class MissingFileException(GRAFIMOException):
    """
        Raise an exception when a file is missing
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        repr(self.message)

class ValueException(GRAFIMOException):
    """
        Raise an exception when a wrong value type is given to
        a function
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        repr(self.message)

class ScaledScoreMatrixException(GRAFIMOException):
    """
        Raise an exception if there are errors regarding the scaled
        score matrix of the motif
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        repr(self.message)

class WrongPathException(GRAFIMOException):
    """
        Raise an exception if is given a wrong path
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        repr(self.message)

class SubprocessException(GRAFIMOException):
    """
        Raise an exception if a subprocess call returns exit status
        different from 0
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        repr(self.message)

class VGException(GRAFIMOException):
    """
        Raise an exception for errors regarding vg tool
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        repr(self.message)

class PipelineException(GRAFIMOException):
    """
        Raise an exception when an error is encountered during the
        pipeline choice or initialization
    """

    def __init__(self, message):

        self.message = message

    def __str__(self):

        repr(self.message)

