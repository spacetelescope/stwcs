import os

from stsci.tools import teal
from .. import __version__
from ..wcsutil import headerlet

__taskname__ = __name__.split('.')[-1] # needed for help string
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

    if configObj['output'] in ['', ' ', None]:
        print('=' * 60)
        print('ERROR:')
        print('    No valid "output" parameter value provided!')
        print('    Please restart this task and provide a value for this parameter.')
        print('=' * 60)
        return

    # create dictionary of remaining parameters, deleting extraneous ones
    # such as those above
    cdict = configObj.dict()
    # remove any rules defined for the TEAL interface
    if "_RULES_" in cdict: del cdict['_RULES_']
    del cdict['_task_name_']
    del cdict['filename']
    del cdict['output']

    # Call function with properly interpreted input parameters
    # Syntax: extract_headerlet(filename, output, extnum=None, hdrname=None,
    #                clobber=False, verbose=100)
    headerlet.extract_headerlet(configObj['filename'], configObj['output'],
                                **cdict)
