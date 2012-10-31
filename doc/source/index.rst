.. STWCS documentation master file, created by
   sphinx-quickstart on Thu Sep 23 09:18:47 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to STWCS's documentation!
=================================

This package provides support for WCS based distortion models and coordinate 
transformation. It relies on PyWCS (based on WCSLIB). It consists of two subpackages:
UPDATEWCS and WCSUTIL. UPDATEWCS performs corrections to the basic WCS and includes 
other distortion infomation in the science files as header keywords or file extensions.
WCSUTIL provides an HSTWCS object which extends pywcs.WCS object and provides HST instrument
specific information as well as methods for coordinate tarnsformaiton. WCSUTIL also provides 
functions for manipulating alternate WCS descriptions in the headers.

Contents:

.. toctree::
   :maxdepth: 2
   
   wcsutil
   updatewcs
   headerlet
      
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

