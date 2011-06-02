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

__docformat__ = 'restructuredtext'

DEGTORAD = fileutil.DEGTORAD
RADTODEG = fileutil.RADTODEG

__version__ = '0.8'

try:
    import svn_version
    __svn_version__ = svn_version.__svn_version__
except ImportError:
    __svn_version__ = 'Unable to determine SVN revision'