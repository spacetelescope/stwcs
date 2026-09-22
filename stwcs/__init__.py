"""STWCS

This package provides support for WCS based distortion models and coordinate
transformation. It relies on astropy.wcs (based on WCSLIB). It consists of
two subpackages:

* updatewcs: Performs corrections to the basic WCS and includes other distortion
  infomation in the science files as header keywords or file extensions.
* wcsutil:  Provides an HSTWCS object which extends astropy.wcs.WCS object and
  provides HST instrument specific information as well as methods for coordinate
  transformation. wcsutil also provides functions for manipulating alternate WCS
  descriptions in the headers.

"""
import logging

from . import distortion  # noqa
from .version import __version__

# Root logger for the whole package. All submodules use
# ``logging.getLogger(__name__)`` so their loggers (e.g. "stwcs.wcsutil.altwcs")
# are children of this "stwcs" logger and propagate up to it. A NullHandler is
# attached here (standard library practice) so nothing is emitted unless the
# application configures logging itself, e.g.:
#
#     import logging
#     handler = logging.StreamHandler()
#     handler.setFormatter(logging.Formatter("%(name)s - %(levelname)s - %(message)s"))
#     logging.getLogger("stwcs").addHandler(handler)
#     logging.getLogger("stwcs").setLevel(logging.DEBUG)
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())
