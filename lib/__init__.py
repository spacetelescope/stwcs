""" STWCS

This package provides support for WCS based distortion models and coordinate 
transformation. It relies on PyWCS (based on WCSLIB). It consists of two subpackages:
updatewcs and wcsutil. Updatewcs performs corrections to the basic WCS and includes 
other distortion infomation in the science files as header keywords or file extensions.
Wcsutil provides an HSTWCS object which extends pywcs.WCS object and provides HST instrument
specific information as well as methods for coordinate tarnsformaiton. Wcsutil also provides 
functions for manipulating alternate WCS descriptions in the headers.

"""
from __future__ import division # confidence high

import distortion
import pywcs
from pytools import fileutil

import logging
log_filename = 'stwcs.log'
logging.basicConfig(filename=log_filename, level=logging.DEBUG, filemode='w')
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger('STWCS')
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(log_filename)
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)
__docformat__ = 'restructuredtext'

DEGTORAD = fileutil.DEGTORAD
RADTODEG = fileutil.RADTODEG

__version__ = '0.8'

try:
    import svn_version
    __svn_version__ = svn_version.__svn_version__
except ImportError:
    __svn_version__ = 'Unable to determine SVN revision'