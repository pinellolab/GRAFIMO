"""Custom Exceptions for GRAFIMO and exception hook handler definition.

The base class GRAFIMOException inherits from Exception. 
GRAFIMOException is extended with custom exceptions, handling errors
occurring during GRAFIMO run.
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


class FileReadError(GRAFIMOException):
    """Raise when errors occurred during generic file reading.
    """
    pass


class FileWriteError(GRAFIMOException):
    """Raise when errors occurred during generic file writing.
    """
    pass


class FileFormatError(GRAFIMOException):
    """Raise when errors occurred during generic file format checks.
    """
    pass


class VGError(GRAFIMOException):
    """Raise when errors occurred during calls to VG
    """
    pass


class MotifFileFormatError(GRAFIMOException):
    """Raise when motif PWM file format is not recognized by GRAFIMO.
    """
    pass


class MotifFileReadError(GRAFIMOException):
    """Raise when errors occurred during motif PWM file parsing.
    """
    pass


class BGFileError(GRAFIMOException):
    """Raise when errors occurred during background distribution file parsing or
    during background distribution computation.
    """
    pass


class MotifProcessingError(GRAFIMOException):
    """Raise when errors occurred during motif position weight matrix processing
    steps.
    """ 
    pass


# ------------------------------------------------------------------------------
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





class NotValidAlphabetException(GRAFIMOException):
    """Raise when a not valid DNA alphabet is given
    """

    pass


class NotValidFFException(GRAFIMOException):
    """Raise if the current file has a not allowed file format
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


class PipelineException(GRAFIMOException):
    """Raise when an errors are encountered during the pipeline choice 
    or initialization
    """

    pass

