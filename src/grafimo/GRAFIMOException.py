"""Custom Exception redefinition for GRAFIMO.

The base class GRAFIMOException, which inherits from Exception, is 
then extended with custome exception fitting possible errors which can
occurr while running GRAFIMO analysis

"""

class GRAFIMOException(Exception):
    """Class to represent an excpetion which can occur ehile running 
    GRAFIMO.

    It is the base for the other more specific GRAFIMO exceptions.
    """

    pass


class DependencyError(GRAFIMOException):
    """Raise if one of the external dependendencies needed by GRAFIMO 
    cannot be found.
    """

    pass


class NoDataFrameException(GRAFIMOException):
    """Raise when a given object is not an instance of pandas.DataFrame
    type
    """

    pass


class WrongMotifWidthException(GRAFIMOException):
    """Raise when the motif is initialized with a not valid width, such 
    as negative or floating point measure
    """

    pass


class WrongMotifIDException(GRAFIMOException):
    """Raise when the motif object is initialized with a wrong motif 
    ID
    """

    pass


class WrongMotifNameException(GRAFIMOException):
    """Raise when the motif object is initialized with a wrong motif 
    name
    """

    pass


class NotValidMotifMatrixException(GRAFIMOException):
    """Raise when the motif object is initialized with a not valid motif 
    matrix
    """

    pass


class NotValidBGException(GRAFIMOException):
    """Raise when a not valid background distribution is given to the 
    motif object
    """

    pass


class NotValidAlphabetException(GRAFIMOException):
    """Raise when a not valid DNA alphabet is given
    """

    pass


class NotValidFFException(GRAFIMOException):
    """Raise if the current file has a not allowed file format
    """

    pass


class FileReadingException(Exception):
    """Raise if is not possible to read the current file
    """

    pass


class ValueException(GRAFIMOException):
    """Raise when a wrong value type is given
    """
    
    pass


class ScaledScoreMatrixException(GRAFIMOException):
    """Raise errors happens computing or loading the scaled score matrix 
    of the motif
    """

    pass


class WrongPathException(GRAFIMOException):
    """Raise if is given a wrong path
    """

    pass


class SubprocessError(GRAFIMOException):
    """Raise if a subprocess call returned an exit status different 
    from 0 (an error occurred)
    """

    pass


class VGException(GRAFIMOException):
    """Raise when errors occur while calling VG
    """

    pass


class PipelineException(GRAFIMOException):
    """Raise when an errors are encountered during the pipeline choice 
    or initialization
    """

    pass

