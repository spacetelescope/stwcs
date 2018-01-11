.. _updatewcs:

**************************************
updatewcs
**************************************
UPDATEWCS applies corrections to the WCS of an HST science file
and adds reference information as header keywords and fits
file extensions so that a science file contains all necessary
information to represent astrometrically precise positions.
The order in which the corrections are applied is important
and is as follows:

- Detector to Image Correction

- Apply Time dependent distortion (if applciable)

- Recomputing the basic WCS

- Apply Velocity Aberration Correction

- Apply polynomial distortion through the SIP coefficients

- Apply non-polynomial distortion

Mathematically the entire transformation from detector to sky
coordinates is described by:

.. math:: (x',y') &= DET2IM(x,y)

.. math:: \binom{u'}{v'} &= \binom{x' - CRPIX1}{y' - CRPIX2}

.. math:: \left( \begin{array}{ll}
         \alpha \\
         \delta \\
         \end{array} \right) &=
      \left( \begin{array}{ll}
      CRVAL1 \\
      CRVAL2\\
      \end{array} \right) +
      \left( \begin{array}{cc}
      CD11 & CD12 \\
      CD21 & CD22\\
      \end{array} \right)
      \left( \begin{array}{ll}
      u' + f(u',v') + LT_x(x',y') \\
      v' + g(u',v') + LT_y(x',y') \\
      \end{array} \right)

where f(u',v') and g(u',v') represent the polynomial distortion correction specified as

.. math:: f(u',v') = \sum_{p+q=2}^{AORDER} A_{pq} {u'}^{p} {v'}^{q}
          \\
          g(u',v')  = \sum_{p+q=2}^{BORDER} B_{pq} {u'}^{p} {v'}^{q}


where

* x', y' are the initial coordinates x,y with the 68th column correction applied through the DET2IM convention
* u',v' are the DET2IM-corrected coordinates relative to CRPIX1,CRPIX2
* LT<sub>x</sub>, LT<sub>y</sub> is the residual distortion in the lookup tables written to the header using the FITS Distortion Paper lookup table convention
* A, B are the SIP coefficients specified using the SIP convention

.. toctree::
   :maxdepth: 2

   updatewcs_ui
   wcs_corrections
   utils
   astrometry_utils

