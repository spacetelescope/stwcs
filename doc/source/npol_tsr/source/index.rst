.. NPOL Reference File TSR master filee.

TSR 2012-03: NPOL Reference File
================================

   Abstract:
   Authors: Nadezhda Dencheva, Warren Hack
   Date: 12 Oct 2010

   The HST pipeline originally used two types of reference files to correct for distortion:
   ``IDCTAB`` files which contain the coefficients for a polynomial correction and
   ``DGEO`` files which are images with the residual distortion. A new format for the residual
   distortion, called ``NPOL`` files (extension ``_npl.fits``),  is presented in this document.
   The conversion from ``DGEO`` files to ``NPOL`` files is described and an example of the format is given
   using ACS/WFC F606W filter. Tests of the conversion procedure show that the
   differences between the ``DGEO`` files and ``NPOL`` files are of the order of :math:`10^{-4}`
   pixels.

Contents:

.. toctree::
   :maxdepth: 2

   npol_tsr
