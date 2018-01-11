Merging Of The Conventions
==========================

The full implementation of all these elements ends up merging the SIP, DET2IM and FITS Distortion Paper 
conventions to create a new version of the figure from the FITS Distortion Paper which illustrates the conversion
of detector coordinates to world coordinates. This implementation works in the following way: 

 #. Apply detector to image correction (DET2IM) to input pixel values
 #. Apply SIP coefficients to DET2IM-corrected pixel values
 #. Apply lookup table correction to DET2IM-corrected pixel values
 #. Add the results of the SIP and lookup table corrections
 #. Apply the WCS transformation in the CD matrix to the summed results to get the intermediate world coordinates
 #. Add the CRVAL keyword values to the transformed positions to get the final world coordinates 
   
The computations to perform these steps can be described approximately using: 

.. math:: (x',y') &= DET2IM(x,y)
   :label: Equation 2

.. math:: \binom{u'}{v'} &= \binom{x' - CRPIX1}{y' - CRPIX2}
   :label: Equation 3

.. _equation4:

.. math:: 
   :label: Equation 4
   
      \left( \begin{array}{ll}
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

.. math:: 
   :label: Equation 5

    f(u',v') = \sum_{p+q=2}^{AORDER} A_{pq} {u'}^{p} {v'}^{q} 
    \\
    g(u',v')  = \sum_{p+q=2}^{BORDER} B_{pq} {u'}^{p} {v'}^{q}


where

* x', y' are the initial coordinates x,y with the 68th column correction applied 
  through the DET2IM convention
* u',v' are the DET2IM-corrected coordinates relative to CRPIX1,CRPIX2
* :math:`LT_{x}, LT_{y}` is the residual distortion in the lookup tables 
  written to the header using the FITS Distortion Paper lookup table convention
* A, B are the SIP coefficients specified using the SIP convention

These equations do not take into account the deprojection from the tangent plane to 
sky coordinates. The complete Detector To Sky Coordinate Transformation is based on 
the CTYPE keyword. 

.. _figure3:

.. figure:: /images/pipeline.png

   Coordinate Transformation Pipeline

