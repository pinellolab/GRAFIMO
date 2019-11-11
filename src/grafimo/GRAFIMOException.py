"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

The script contains definitions of custom Exceptions
for GRAFIMO

"""

class NoDataFrameException(Exception):
    """
        Raise when a given object is not an instance of pandas.DataFrame
    """
    def __init__(self, message):

        self.message=message

    def __str__(self):

        return repr(self.message)

class WrongMotifWidthException(Exception):
    """
        Raise an exception when the motif is initialized with a
        not valid width
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        return repr(self.message)

class WrongMotifIDException(Exception):
    """
        Raise an exception when the motif object is initialized
        with a wrong ID
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        return repr(self.message)

class WrongMotifNameException(Exception):
    """
        Raise an exception when the motif object is initialized
        with a wrong name
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        return repr(self.message)

class NotValidMotifMatrixException(Exception):
    """
        Raise an exception when is given a not valid motif matrix
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        return repr(self.message)

class NotValidBGException(Exception):
    """
        Raise  an exception when a not valid background distribution is
        given to the motif object
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        return repr(self.message)

class NotValidAlphabetException(Exception):
    """
        Raise an exception when a not valid alphabet is given
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        return repr(self.message)

class NotValidFFException(Exception):
    """
        Raise an exception if it is given a file in a not allowed
        file format
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        return repr(self.message)

class FileReadingException(Exception):
    """
        Raise an exception if something went wrong during a file
        reading
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        return repr(self.message)

class MissingFileException(Exception):
    """
        Raise an exception when a file is missing
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        repr(self.message)

class ValueException(Exception):
    """
        Raise an exception when a wrong value type is given to
        a function
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        repr(self.message)

class ScaledScoreMatrixException(Exception):
    """
        Raise an exception if there are errors regarding the scaled
        score matrix of the motif
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        repr(self.message)

class WrongPathException(Exception):
    """
        Raise an exception if is given a wrong path
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        repr(self.message)

class SubprocessException(Exception):
    """
        Raise an exception if a subprocess call returns exit status
        different from 0
    """

    def __init__(self, message):

        self.message=message

    def __str__(self):

        repr(self.message)