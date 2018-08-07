import os

from astropy.io import fits
from stsci.tools import parseinput
from stsci.tools import fileutil
from stsci.tools import teal
from .. import __version__
from .. import updatewcs


allowed_corr_dict = {'vacorr': 'VACorr', 'tddcorr': 'TDDCorr', 'npolcorr': 'NPOLCorr',
                     'd2imcorr': 'DET2IMCorr'}


__taskname__ = __name__.split('.')[-1] # needed for help string
__package__ = updatewcs.__name__
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
            helpString += eval('.'.join([__package__, __taskname__, '__doc__']))

    else:
        helpString = 'file://' + htmlfile

    return helpString


def run(configObj=None):

    # Interpret primary parameters from configObj instance
    input = configObj['input']

    # create dictionary of remaining parameters, deleting extraneous ones
    # such as those above
    cdict = configObj.dict()
    # remove any rules defined for the TEAL interface
    if "_RULES_" in cdict: del cdict['_RULES_']
    del cdict['_task_name_']
    del cdict['input']
    del cdict['extname']

    # parse input
    input, altfiles = parseinput.parseinput(configObj['input'])

    # Insure that all input files have a correctly archived
    #    set of OPUS WCS keywords
    # Legacy files from OTFR, like all WFPC2 data from OTFR, will only
    #   have the OPUS WCS keywords archived using a prefix of 'O'
    # These keywords need to be converted to the Paper I alternate WCS
    #   standard using a wcskey (suffix) of 'O'
    # If an alternate WCS with wcskey='O' already exists, this will copy
    #   the values from the old prefix-'O' WCS keywords to insure the correct
    #   OPUS keyword values get archived for use with updatewcs.
    #
    for file in input:
        # Check to insure that there is a valid reference file to be used
        idctab = fits.getval(file, 'idctab')
        if not os.path.exists(fileutil.osfn(idctab)):
            print('No valid distortion reference file ', idctab, ' found in ', file, '!')
            raise ValueError

    # Re-define 'cdict' to only have switches for steps supported by that instrument
    # the set of supported steps are defined by the dictionary
    #    updatewcs.apply_corrections.allowed_corrections
    #
    for file in input:
        # get instrument name from input file
        instr = fits.getval(file, 'INSTRUME')
        # make copy of input parameters dict for this file
        fdict = cdict.copy()
        # Remove any parameter that is not part of this instrument's allowed corrections
        for step in allowed_corr_dict:
            if allowed_corr_dict[step] not in updatewcs.apply_corrections.allowed_corrections[instr]:
                fdict[step]
        # Call 'updatewcs' on correctly archived file
        updatewcs.updatewcs(file, **fdict)
