import os
from stsci.tools import teal
from stsci.tools import parseinput

from .. import __version__
from ..wcsutil import headerlet

__taskname__ = __name__.split('.')[-1]  # needed for help string
__package__ = headerlet.__name__
# __version__ = stwcs.__version__
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

    if configObj['hdrname'] == '' and configObj['hdrext'] is None and \
            configObj['distname'] == '':
        print('=' * 60)
        print('ERROR:')
        print('    No valid "hdrname", "hdrext" or "distname" parameter value provided!')
        print('    Please restart this task and provide a value for one of these parameters.')
        print('=' * 60)
        return
    filename = parseinput.parseinput(configObj['filename'])[0]
    # Call function with properly interpreted input parameters
    # Syntax: delete_headerlet(filename, hdrname=None, hdrext=None, distname=None)
    headerlet.delete_headerlet(filename,
                               hdrname=configObj['hdrname'],
                               hdrext=configObj['hdrext'],
                               distname=configObj['distname'],
                               logging=configObj['logging'])
