.. STWCS documentation master file, created by
   sphinx-quickstart on Thu Sep 23 09:18:47 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

STWCS documentation
===================

This package provides support for WCS based distortion models and coordinate
transformation developed to support HST imaging observations.  It consists of two subpackages -
`updatewcs <updatewcs.html>`_ and `wcsutil <wcsutil.html>`_. `updatewcs  <updatewcs.html>`_
performs corrections to the basic WCS and inserts
other distortion infomation in the science files as header keywords or file extensions.
`wcsutil <wcsutil.html>`_ provides an `wcsutil.STWCS` object which extends
astropy's `astropy.wcs.WCS` provides HST
instrument specific information as well as methods for coordinate transformation. `wcsutil <wcsutil.html>`_ also provides
functions for manipulating alternate WCS descriptions in the headers.
`STWCS <https://github.com/spacetelescope/stwcs>`__ is based on
`astropy.wcs <http://docs.astropy.org/en/stable/wcs/index.html>`_.


Contents:

.. toctree::
   :maxdepth: 2

   wcsutil
   updatewcs
   headerlet
   technical_reports

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

