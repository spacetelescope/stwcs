from __future__ import division # confidence high

#import all needed modules here to avoid relative imports
#import mappings
import utils
import distortion
import pywcs
from pytools import fileutil

DEGTORAD = fileutil.DEGTORAD
RADTODEG = fileutil.RADTODEG

try:
    import svn_version
    __svn_version__ = svn_version.__svn_version__
except ImportError:
    __svn_version__ = 'Unable to determine SVN revision'