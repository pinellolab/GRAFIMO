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


class DependencyException(GRAFIMOException):
    """
        raise an excepton if one of the external
        dependendencies needed by GRAFIMO is
        not found
    """

    pass


class NoDataFrameException(GRAFIMOException):
    """
        Raise when a given object is not an instance of pandas.DataFrame
    """

    pass


class WrongMotifWidthException(GRAFIMOException):
    """
        Raise an exception when the motif is initialized with a
        not valid width
    """

    pass


class WrongMotifIDException(GRAFIMOException):
    """
        Raise an exception when the motif object is initialized
        with a wrong ID
    """

    pass


class WrongMotifNameException(GRAFIMOException):
    """
        Raise an exception when the motif object is initialized
        with a wrong name
    """

    pass


class NotValidMotifMatrixException(GRAFIMOException):
    """
        Raise an exception when is given a not valid motif matrix
    """

    pass


class NotValidBGException(GRAFIMOException):
    """
        Raise  an exception when a not valid background distribution is
        given to the motif object
    """

    pass


class NotValidAlphabetException(GRAFIMOException):
    """
        Raise an exception when a not valid alphabet is given
    """

    pass


class NotValidFFException(GRAFIMOException):
    """
        Raise an exception if it is given a file in a not allowed
        file format
    """

    pass


class FileReadingException(Exception):
    """
        Raise an exception if something went wrong during a file
        reading
    """

    pass


class MissingFileException(GRAFIMOException):
    """
        Raise an exception when a file is missing
    """

    pass


class ValueException(GRAFIMOException):
    """
        Raise an exception when a wrong value type is given to
        a function
    """
    
    pass


class ScaledScoreMatrixException(GRAFIMOException):
    """
        Raise an exception if there are errors regarding the scaled
        score matrix of the motif
    """

    pass


class WrongPathException(GRAFIMOException):
    """
        Raise an exception if is given a wrong path
    """

    pass


class SubprocessException(GRAFIMOException):
    """
        Raise an exception if a subprocess call returns exit status
        different from 0
    """

    pass


class VGException(GRAFIMOException):
    """
        Raise an exception for errors regarding vg tool
    """

    pass


class PipelineException(GRAFIMOException):
    """
        Raise an exception when an error is encountered during the
        pipeline choice or initialization
    """

    pass
