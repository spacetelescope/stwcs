SIP Convention
==============

Current implementations of distortion models in FITS headers have been limited to simply 
describing polynomial models. The prime example of this would be the implementation of SIP 
in WCSTOOLS and recognition of SIP keywords by DS9 as used for Spitzer data [SIPConvention]_. 
The new keywords defined by the SIP standard and used by PyWCS are::

 A_ORDER =   n  / polynomial order, axis 1, detector to sky
 A_i_j          / High order coefficients for X axis
 B_ORDER =   m  /  polynomial order, axis 2, detector to sky
 B_i_j          / High order coefficients for axis 2

These SIP keywords get used in conjunction with the linear WCS keywords defined
with these values::

 CTYPE1  = 'RA---TAN-SIP'
 CTYPE2  = 'DEC--TAN-SIP'
 CDi_j          / Linear terms of distortion plus scale and orientation

The SIP convention retains the use of the current definition of the CD matrix where the 
linear terms of the distortion model are folded in with the orientation and scale at the 
reference point for each chip to provide the best linear approximation to the distortion 
available. The SIP convention gets applied to the input pixel positions by applying the 
higher-order coefficients A_i_j, B_i_j, then by applying the CD matrix and adding the CRVAL 
position to get the final world coordinates.

This convention was created from the original form of the FITS Distortion Paper standards, but the 
FITS Distortion Paper proposal has since changed to use a different set of keywords and conventions. 

A sample ACS/WFC SCI header can be found in :ref:`appendix1` to illustrate how these 
keywords actually get populated for an image.  The current implementation does not 
take advantage of the A_DMAX, B_DMAX, SIPREFi or SIPSCLi keywords, so these keywords
are not written out to the SCI header.

Velocity Aberration Correction
------------------------------

This correction simply serves as a correction to the overall linear scale of the field of view
due to velocity aberration observed due to the motion of HST in orbit.  The typical plate scale
for HST cameras results in a measurable velocity aberration with variations from the center of
the field of view to the edge on the order of 0.1 pixels. More details about this correction can
be found in `Appendix A.3 of the DrizzlePac Handbook
<https://hst-docs.stsci.edu/drizzpac>`_.

This scale factor gets computed by
the HST ground systems for start of each exposure and recorded as the VAFACTOR keyword in each
image's science extension header. This term, though, does not get included in the default 
CD matrix computed by the ground systems. As a result, it needs to be accounted for when reading in the 
distortion model polynomial coefficients from the IDCTAB reference table. The VAFACTOR scaling factor
gets folded into the computation of new values for the CD matrix for this specific exposure 
without requiring any further use of the VAFACTOR keyword when applying this distortion 
model to the science image. It also gets used to correct the reference position of each chip
on the sky, each chip's CRVAL value, to account for this aberration. 


Time-Dependent Distortion
-------------------------

Calibration of HST/ACS imaging data required the addition of a time dependent skew in addition
to the other distortion terms.  This skew represented a linear correction to the polynomial model
and its residuals.  This correction gets applied to the polynomial coefficients and
the residuals from the polynomial model when they are evaluated for each image.  As a result, the 
SIP keywords as written out to each HST/ACS image header
reflects this time-dependent correction without the need for any further evaluation of this skew.

