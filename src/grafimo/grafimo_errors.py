"""Custom Exceptions for GRAFIMO and exception hook handler definition.

The base class GRAFIMOError inherits from Exception. 
GRAFIMOError is extended with custom exceptions, handling errors
occurring during GRAFIMO run.
"""


class GRAFIMOError(Exception):
    """Class to represent an excpetion which can occur ehile running 
    GRAFIMO.

    It is the base for the other more specific GRAFIMO exceptions.
    """
    pass


class DependencyError(GRAFIMOError):
    """Raise if one of the external dependendencies needed by GRAFIMO 
    cannot be found.
    """
    pass


class FileReadError(GRAFIMOError):
    """Raise when errors occurred during generic file reading.
    """
    pass


class FileWriteError(GRAFIMOError):
    """Raise when errors occurred during generic file writing.
    """
    pass


class FileFormatError(GRAFIMOError):
    """Raise when errors occurred during generic file format checks.
    """
    pass


class VGError(GRAFIMOError):
    """Raise when errors occurred during calls to VG
    """
    pass


class MotifFileFormatError(GRAFIMOError):
    """Raise when motif PWM file format is not recognized by GRAFIMO.
    """
    pass


class MotifFileReadError(GRAFIMOError):
    """Raise when errors occurred during motif PWM file parsing.
    """
    pass


class BGFileError(GRAFIMOError):
    """Raise when errors occurred during background distribution file parsing or
    during background distribution computation.
    """
    pass


class MotifProcessingError(GRAFIMOError):
    """Raise when errors occurred during motif position weight matrix processing
    steps.
    """ 
    pass


class NotValidMotifMatrixError(GRAFIMOError):
    """Raise when the motif object is initialized with a not valid motif 
    matrix
    """
    pass


class SubprocessError(GRAFIMOError):
    """Raise if a subprocess call returned an exit status different 
    from 0 (an error occurred)
    """
    pass

