import os

from stsci.tools import teal
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
    Return useful help from a file in the script directory called __taskname__.help
    """
    install_dir = os.path.dirname(__file__)
    htmlfile = os.path.join(install_dir, 'htmlhelp', __taskname__ + '.html')
    helpfile = os.path.join(install_dir, __taskname__ + '.help')
    if docstring or (not docstring and not os.path.exists(htmlfile)):
        helpString = __taskname__ + ' Version ' + __version__ + '\n\n'
        if os.path.exists(helpfile):
            helpString += teal.getHelpFileAsString(__taskname__, __file__)
        else:
            helpString += headerlet.archive_as_headerlet.__doc__

    else:
        helpString = 'file://' + htmlfile

    return helpString


def run(configObj=None):

    if configObj['hdrname'] in ['', ' ', 'INDEF']:
        print('=' * 60)
        print('ERROR:')
        print('    No valid "hdrname" parameter value provided!')
        print('    Please restart this task and provide a value for this parameter.')
        print('=' * 60)
        return

    str_kw = ['wcsname', 'destim', 'sipname', 'npolfile', 'd2imfile',
              'descrip', 'history', 'author']

    # create dictionary of remaining parameters, deleting extraneous ones
    # such as those above
    cdict = configObj.dict()
    # remove any rules defined for the TEAL interface
    if "_RULES_" in cdict: del cdict['_RULES_']
    del cdict['_task_name_']
    del cdict['filename']
    del cdict['hdrname']

    # Convert blank string input as None
    for kw in str_kw:
        if cdict[kw] == '': cdict[kw] = None
    if cdict['wcskey'].lower() == 'primary': cdict['wcskey'] = ' '

    # Call function with properly interpreted input parameters
    # Syntax: archive_as_headerlet(filename, sciext='SCI', wcsname=None, wcskey=None,
    #                    hdrname=None, destim=None,
    #                    sipname=None, npolfile=None, d2imfile=None,
    #                    author=None, descrip=None, history=None,
    #                    hdrlet=None, clobber=False)
    headerlet.archive_as_headerlet(configObj['filename'], configObj['hdrname'],
                                   **cdict)
