import os

from stsci.tools import teal
from ..wcsutil import headerlet
import stwcs

__taskname__ = __name__.split('.')[-1]  # needed for help string
__package__ = headerlet.__name__
__version__ = stwcs.__version__
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
            helpString += eval('.'.join([__package__, __taskname__, '__doc__']))
    else:
        helpString = 'file://' + htmlfile

    return helpString


def run(configObj=None):

    headerlet.attach_headerlet(configObj['filename'], configObj['hdrlet'],
                               configObj['logging'])
