"""GRAFIMO redefinition of ArgumentParser and HelpFormatter

The ArgumentParser and HelpFormatter classes are customized to 
assign a different style to the help with respect to the default
one and when a wrong argument is given GRAFIMO reminds the user
how to call the command-line help. 
"""

from argparse import ArgumentParser, SUPPRESS, HelpFormatter
from grafimo.grafimo import __version__
from typing import List, Tuple, Dict, Optional
import sys


class GRAFIMOArgumentParser(ArgumentParser):
    """
    This class redefines the ArgumentParser class from argparse module.

    The redefinition is required to change the argparse default definition
    of the help. So, when raising an error by argument parser will be 
    shown the new GRAFIMO message

    ...


    Methods
    -------
    error(msg)
        Prints a new message instead of argparse deafault one when an error
        is raised parsing the arguments

    """

    class GRAFIMOHelpFormatter(HelpFormatter):
        """
        This class defines a new help format for GRAFIMO

        ...


        Methods
        -------
        add_usage(usage, actions, groups, prefix=None)
            style for GRAFIMO help
        """
        #---------------------------------------------------------------
        # GRAFIMOHelpFormatter methods
        #---------------------------------------------------------------
        def add_usage(
            self, 
            usage: str, 
            actions: str, 
            groups: str, 
            prefix: Optional[str] = 'None'
        ) -> None:
        
            if usage is not SUPPRESS:
                args = usage, actions, groups, ''
                self._add_item(self._format_usage, args)

    # end of GRAFIMOHelpFormatter

    #-------------------------------------------------------------------
    # GRAFIMOArgumentParser methods
    #-------------------------------------------------------------------
    def __init__(
        self, 
        *args: Tuple, 
        **kwargs: Dict
    ):
        """
        Parameters
        ----------
        args : tuple
            Arguments to initialize the ArgumentParser object
        kwargs : dict
            Help format restyling arguments
        """

        kwargs['formatter_class'] = self.GRAFIMOHelpFormatter
        kwargs['usage'] = kwargs['usage'].replace("{version}", __version__)
        super().__init__(*args, **kwargs)


    def error(
        self, 
        msg: str
    ) -> None:
        """Print the given message when argparse raise an error.

        At the bottom is also printed a new message telling the 
        user how to call the halp

        Parameters
        ----------
        msg : str
            Message which will be shown when raising an error
        """
        
        sys.stderr.write("Run 'grafimo --help' to see the usage\n")
        self.exit(2, "\n{0}: ERROR: {1}\n".format(self.prog, msg))

# end of GRAFIMOArgumentParser

