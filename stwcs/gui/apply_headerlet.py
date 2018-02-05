import os
from stsci.tools import teal, parseinput

from .. import __version__
from ..wcsutil import headerlet

__taskname__ = __name__.split('.')[-1]  # needed for help string
__package__ = headerlet.__name__
# __version__ = stwcs.__version__


############### Interfaces used by TEAL ###############

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

    # start by interpreting filename and hdrlet inputs
    filename = parseinput.parseinput(configObj['filename'])[0]
    hdrlet = parseinput.parseinput(configObj['hdrlet'])[0]

    if configObj['primary']:
        # Call function with properly interpreted input parameters
        # Syntax: apply_headerlet_as_primary(filename, hdrlet, attach=True,
        #            archive=True, force=False, verbose=False)
        headerlet.apply_headerlet_as_primary(filename,
                                             hdrlet, attach=configObj['attach'],
                                             archive=configObj['archive'], force=configObj['force'],
                                             logging=configObj['logging'])
    else:
        wcsname = configObj['wcsname']
        if wcsname in ['', ' ', 'INDEF']: wcsname = None
        wcskey = configObj['wcskey']
        if wcskey == '': wcskey = None
        # Call function with properly interpreted input parameters
        #         apply_headerlet_as_alternate(filename, hdrlet, attach=True,
        #                        wcskey=None, wcsname=None, verbose=False)
        headerlet.apply_headerlet_as_alternate(filename,
                                               hdrlet, attach=configObj['attach'],
                                               wcsname=wcsname, wcskey=wcskey,
                                               logging=configObj['logging'])
