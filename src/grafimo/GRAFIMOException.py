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


class MotifMatrixError(GRAFIMOException):
    """Raise when invalid motif matrix is used
    """
    pass


class SubprocessError(GRAFIMOException):
    """Raise if a subprocess call returned an exit status different 
    from 0 (an error occurred)
    """
    pass
