.. _astrometry_utils:

**************************************
Astrometry Database Interface
**************************************
This module provides an interface for accessing a database which contains
updated WCS solutions for many images in the HST archive.  The exact set of
images with new solutions in this database will change over time, with the initial
release containing updates only for ACS and WFC3 direct images.  This database
will be designed to provide multiple WCS solutions for each exposure depending
on what astrometric frame was used to define the new WCS solution; for example,
an image may be matched to Pan-STARRS sources in one solution and to GAIA sources
in another solution.

Each new WCS solution will be appended to the exposure by this module as a
headerlet extension, taking advantage of the headerlet definition provided by the
STWCS package, found in :ref:`headerlet` documentation.
The headerlet library allows a user
to manage all appended solutions and determine which one gets used by the
exposure at anytime, since only 1 WCS can be used at a time.

Usage
------
This interface defines an object that provides the services for updating an
exposure with new WCS solutions. Only 2 calls are typically necessary to update
any exposure or set of exposures.

0. Import the package
::

  >>> from stwcs.updatewcs import astrometry_utils as stwcsau

1. Establish a working connection with the database
::

  >>> db = stwcsau.AstrometryDB()

2. For each exposure, query the database and update with any new WCS solutions
::

  >>> db.updateObs(filename)

The environment variables used by this module can be defined by the user to change
what database is being accessed, whether to raise exceptions instead of warnings
when problems arise in connecting to the database or otherwise updating an
exposure, or to completely turn on/off use of this module.  These allow external
code, such as updatewcs (see :ref:`updatewcs`), to include this code and
control it as necessary.


API
----
.. automodule:: stwcs.updatewcs.astrometry_utils

.. autoclass:: stwcs.updatewcs.astrometry_utils.AstrometryDB
  :members:
  :inherited-members:
  :show-inheritance:

.. autofunction:: stwcs.updatewcs.astrometry_utils.apply_astrometric_updates
.. autofunction:: stwcs.updatewcs.astrometry_utils.find_gsc_offset

The :py:func:`find_gsc_offset` function uses a web service to
provide an improved position for both the reference guide star and the secondary (roll) guide star
used in an HST observation. The Guide Star Catalog astrometry has improved a lot over the years.
The original version of GSC (GSC-1) had errors in astrometry as large as 1 or 2 arcseconds.
GSC-2 reduced those errors to about 0.25 arcseconds. More recent versions of the catalog
(now up to GSC 3) reduce the errors further using Gaia positions,
typically 0.01 arcsecond or less, depending on the epoch of the observation.
By knowing what the guide star position was thought to be at the time (from the
original catalog), and what the position really is (from the improved catalog),
one can get a more accurate value for where the telescope was actually pointing.

Version
-------

.. autodata:: stwcs.__version__
