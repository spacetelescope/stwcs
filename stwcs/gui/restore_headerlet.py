import os

from stsci.tools import teal

from .. import __version__
from ..wcsutil import headerlet

__taskname__ = __name__.split('.')[-1] # needed for help string
__package__ = headerlet.__name__
#__version__ = stwcs.__version__
#
#### Interfaces used by TEAL
#
def getHelpAsString(docstring=False):
    """
    return useful help from a file in the script directory called __taskname__.help
    """
    install_dir = os.path.dirname(__file__)
    htmlfile = os.path.join(install_dir, 'htmlhelp', __taskname__ + '.html')
    helpfile = os.path.join(install_dir, __taskname__ + '.help')
    if docstring or (not docstring and not os.path.exists(htmlfile)):
        helpString = __taskname__ + ' Version ' + __version__ + '\n\n'
        if os.path.exists(helpfile):
            helpString += teal.getHelpFileAsString(__taskname__, __file__)
    else:
        helpString = 'file://' + htmlfile

    return helpString


def run(configObj=None):

    if configObj['distname'] not in ['', ' ', 'INDEF']:
        # Call function with properly interpreted input parameters
        # Syntax: restore_all_with_distname(filename, distname, primary,
        #           archive=True, sciext='SCI', verbose=False)
        headerlet.restore_all_with_distname(configObj['filename'],
                                            configObj['distname'], configObj['primary'],
                                            archive=configObj['archive'], sciext=configObj['sciext'],
                                            logging=configObj['logging'])
    else:
        # Call function with properly interpreted input parameters
        #         restore_from_headerlet(filename, hdrname=None, hdrext=None,
        #           archive=True, force=False)
        headerlet.restore_from_headerlet(configObj['filename'],
                                         hdrname=configObj['hdrname'], hdrext=configObj['hdrext'],
                                         archive=configObj['archive'], force=configObj['force'],
                                         logging=configObj['logging'])
