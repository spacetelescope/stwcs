from __future__ import division # confidence high

from altwcs import *
from hstwcs import HSTWCS

__docformat__ = 'restructuredtext'


def help():
    print 'How to create an HSTWCS object:\n\n'
    print """ \
    1. Using a pyfits HDUList object and an extension number \n
    Example:\n

    fobj = pyfits.open('some_file.fits')\n
    w = wcsutil.HSTWCS(fobj, 3)\n\n

    2. Create an HSTWCS object using a qualified file name. \n
    Example:\n
    w = wcsutil.HSTWCS('j9irw4b1q_flt.fits[sci,1]')\n\n

    3. Create an HSTWCS object using a file name and an extension number. \n
    Example:\n
    w = wcsutil.HSTWCS('j9irw4b1q_flt.fits', ext=2)\n\n

    4. Create an HSTWCS object from WCS with key 'O'.\n
    Example:\n

    w = wcsutil.HSTWCS('j9irw4b1q_flt.fits', ext=2, wcskey='O')\n\n

    5. Create a template HSTWCS object for a DEFAULT object.\n
    Example:\n
    w = wcsutil.HSTWCS(instrument='DEFAULT')\n\n
    """

