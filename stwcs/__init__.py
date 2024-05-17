""" STWCS

This package provides support for WCS based distortion models and coordinate
transformation. It relies on astropy.wcs (based on WCSLIB). It consists of
two subpackages:  updatewcs and wcsutil.

updatewcs performs corrections to the basic WCS and includes other distortion
infomation in the science files as header keywords or file extensions.

Wcsutil provides an HSTWCS object which extends astropy.wcs.WCS object and
provides HST instrument specific information as well as methods for coordinate
transformation. wcsutil also provides functions for manipulating alternate WCS
descriptions in the headers.

"""
import os

from . import distortion  # noqa
from stsci.tools import fileutil  # noqa
from stsci.tools import teal  # noqa

from .version import __version__

try:
    from . import gui
    teal.print_tasknames(gui.__name__, os.path.dirname(gui.__file__))
    print('\n')
except:
    pass
