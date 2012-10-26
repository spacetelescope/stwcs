===============================================
Distortion Correction in HST FITS Files - DRAFT
===============================================

.. abstract::
   :author: Warren Hack, Nadezhda Dencheva, Andy Fruchter, Perry Greenfield
   :date: 18 Sept 2012
   
   A convention for storing distortion information in HST images was developed 
   and implemented in two software packages - PyWCS and STWCS. These changes 
   allow the development of a WCS based version of Multidrizzle, released
   now as AstroDrizzle, and image alignment software, released as tweakreg. 
   The distribution of WCS solutions is discussed.
 
   
Introduction
============

Calibration of the HST Advanced Camera for Surveys (HST/ACS) distortion requires the use 
of several components to the distortion correction; namely, polynomial coefficients, a correction 
for velocity aberration, a time-dependent skew, a lookup table for non-polynomial terms, and a detector defect 
correction. Each of these terms has been derived as part of the calibration effort to address 
separate aspects of the distortion that affects ACS observations. Ideally, each would be applied
independently in the same manner used for deriving the original calibration reference information, 
with the time-dependent skew being folded into the other terms. However, the software for 
applying the distortion models does not support this option. In fact, there is no clear 
accepted standard for specifying distortion corrections in FITS headers. Instead, there are 
several separate proposals for specifying aspects of the distortion, but none by themselves 
allows us to fully specify the distortion already calibrated for ACS, let alone in a modular, 
efficient manner.

This paper describes a composite implementation of a select set of proposed standards which 
supports all aspects of the distortion models for HST instruments without breaking any of the 
conventions/standards. The rules for merging the proposed standards allow software to be defined 
to apply each aspect of this proposal as a separate option while defining the requirements 
necessary to allow them to work together when specified in the header. As a result, the separate 
components essentially become tools where only those conventions appropriate to the observation 
can be used as needed. 

Problems Introduced by the HST/ACS Distortion 
=============================================

All calibrations for HST observations get recorded and applied through the use of 
reference files, separate files which describe some calibration. The geometric 
distortion typically applied to HST images gets recorded as a polynomial component 
in one reference file, and a pixel-by-pixel correction to the polynomial solution 
in a separate reference file. This method allows the distortion to be corrected to 
an accuracy of better than 0.1 pixels. However, this method requires the user to 
obtain the reference files themselves anytime they want to reprocess the data. The 
size of these reference files (up to 200Mb) makes this an expensive requirement for 
the end user. The alternative would be to include the necessary specification of the 
distortion model in the header of the science image itself, as long as it can be done in a 
manner that does not dramatically increase the size of the FITS file. For reference, 
a typical calibrated ACS/WFC image requires a 168Mb file. Thus, we needed an alternative 
to separate reference files which can be implemented in a very efficient manner within 
the image's FITS headers.

The calibrations also represent separate aspects of the detector and the distortion, 
aspects which logically should remain as separate descriptions in the header. The pixels 
for each CCD do not have the same size across the entire chip. The ACS/WFC CCDs manufacturing 
process resulted in the pixels having different sizes every 68.3 columns, and can be represented 
most efficiently and accurately by a 1-D correction that gets applied to every row in the chip. 
This detector level characterization affects all images readout from the CCD regardless of any 
additional distortion applied to the field of view. Logically, this effect should be kept as a 
separate component of the distortion model that gets applied prior to correcting for any other 
distortion. This represents an example of an effect that is best applied sequentially to the image data.

Additional distortions come as a result of the effect of the optics on the field of view. 
These are generally described by low-order polynomials for most instruments, although for 
ACS, an additional non-polynomial correction needed to be taken into account as well. 
Fortunately, the non-polynomial correction can sub-sampled enough to make it practical 
to include in the image headers directly, a correction of up to a quarter of a pixel in some 
areas of the detector. Both components need to be applied to the data in order to align images 
well enough for subsequent data analysis or cosmic-ray rejection.

These corrections could be combined into a single look-up table, yet it would come at the 
cost of additional errors which may not allow the remaining data analysis or cosmic-ray 
rejection to actually succeed. We also have some instruments where there is only a polynomial 
component, requiring the development of support for a polynomial correction and a look-up 
table anyway.

These requirements on the application of the calibrations to HST data leave us with no 
alternative within current FITS standards. As a result, we developed this set of rules 
which allow us to take advantage of the most appropriate conventions for each separate 
component of the distortion model and combine them in an efficient manner which eliminates 
the need for external reference data.

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

A sample ACS/WFC SCI header can be found in :ref:`Appendix1` to illustrate how these 
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
<http://documents.stsci.edu/hst/HST_overview/documents/DrizzlePac/DrizzlePac.cover.html>`_.

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


FITS Distortion Proposal
=========================

The current FITS Distortion Paper conventions [DistortionPaper]_ provide a mechanism for specifying either a lookup table 
or polynomial model for the distortion of each axis. The standard states in Section 2.1: 

``Note that the prior distortion functions,, operate on pixel coordinates (i.e. p  
rather than p− r ), and that the independent variables of the distortion functions 
are the uncorrected pixel or intermediate pixel coordinates. That is, for example, 
we do not allow the possibility of``

.. math::
   :label: Equation 1

   q'_{3} = q_{3} + \delta_{q_{3}}(q'_{1},q'_{2})

The keywords used for describing these corrections use the syntax given in Table 2 of the FITS Distortion Paper. 
For our purposes, the keywords of interest are those related to lookup tables; namely, 

::

 CPDISja        string    2.4.1 distortion code new Prior distortion function type.
 DPja           record    2.4.2 distortion parameter new Parameter for a prior distortion 
                                  function, for use in an image header
                          
This syntax only provides the option to specify one correction at a time for each 
axis of the image. This precludes being able to use this convention to specify both 
a lookup table and a polynomial model at the same time for the same axis. It does not 
state what should be done if the polynomial has been specified using a different 
convention, for example, the SIP convention. Thus, SIP and FITS Distortion Paper should not be 
seen as mutually exclusive. In fact, they may work together rather naturally since the 
SIP and FITS Distortion Paper conventions both assume the corrections will work on the input pixel 
and add to the output frame. 

The sample header in :ref:`Appendix1` shows how these keywords get populated for
an actual reference file; specifically, an NPOLFILE as described in the next section.


Non-polynomial Residual Correction
==================================

ACS and WFPC2 images used the DGEOFILE reference file to specify the residual
correction in X and Y for each and every pixel in each chip of the observation. These
DGEOFILE reference fiels required up to 168Mb each to cover all chips of each camera 
for ACS/WFC images.  
Distortion residuals have been always been calibrated for ACS by looking at the average 
correction that needs to be applied over each 64x64 pixel section of each chip after applying 
the polynomial coefficients correction. This would normally result in a 64 x 32 array of 
residuals for each 4096 x 2048 chip. These arrays get expanded by one value in each 
dimension to support interpolation all the way to the edge of each chip resulting in 65 x 33 arrays of 
distortion correction data. 

.. _figure1:

.. figure:: /images/npol_vector_text.png
   :width: 95 %
   :alt: ACS/WFC F475W NPOLFILE corrections
   :align: center
   
   This figure illustrates the corrections included in the ACS/WFC F475W non-polynomial
   distortion correction included in the new NPOLFILE reference file. Each vector represents
   the correction for a 64x64 pixel section of each chip.


These look-up tables follow the conventions 
in the WCS FITS Distortion Paper [DistortionPaper]_. 
Record-valued keywords are used to map an image in the science extension 
to a distortion array in the ``WCSDVAR extension``. This new type of FITS keywords has been 
implemented in PyFITS and is fully described in [DistortionPaper]_. Specifically, ``DPj.EXTVER`` in the science 
extension header  maps the science image to the correct ``WCSDVAR`` extension. The dimensionality 
of the distortion array is defined by ``DPj.NAXES``. Keywords ``DPj.AXIS.j`` in the ``SCI`` extension 
header are used for mapping image array axis to distortion array axis. In the keywords above j 
is an integer and denotes the axis number. For example, if distortion array axis 1 corresponds 
to image array axis 1 of  a ``SCI`` extension, then ``DP.1.AXIS.1`` = 1.                           
A full example of the keywords added to a ``SCI`` extension header is presented in :ref:`Appendix1`.

A complete description of the conversion of the DGEOFILE reference data into NPOLFILE reference
files can be found in the report on the :ref:`npolfile-tsr`.


NPOLFILE reference File Format
------------------------------

With the goal of including all distortion reference information directly in the 
science image's FITS file, including the full 168Mb DGEOFILE for ACS/WFC images 
would more than double the size of each input image. A new reference
file based on the sub-sampled calibrations, though, would be small enough to serve as the
basis for a new reference file while also being a more direct use of the calibration
data. This new reference file has been called **NPOLFILE** in the FITS image header, 
so that any original DGEOFILE reference filename can be retained in parallel for 
backwards compatibility with the current software. This reference file also 
has a unique suffix, **_npl.fits**, as another means of identifying it as a new 
reference file separate from the current DGEOFILE files. The header for this new 
reference file also remains very simple, as illustrated in :ref:`Appendix2`.

Applying these corrections starts by reading the two 65 x 33 
arrays into memory with each input ACS/WFC chip WCS (one for 
X offsets and one for Y offsets). Bi-linear interpolation based on the input pixel 
position then gets used on-the-fly to extract the final offset from this reference 
file. Initial versions of these sub-sampled NPOLFILE reference files for ACS have 
been derived from the current full-size DGEOFILEs, and testing indicates residuals 
only on the order of 0.02 pixels or less remain when compared to the original calibration. 

Detector To Image Correction
============================

The last element of the distortion which remains to be described is the fixed column 
(or row) width correction. This needs to be applied as a correction to the input pixel 
position and the output of this correction is to be used as input to the polynomial and 
non-polynomial distortion corrections.

The adopted implementation is based on the FITS Distortion Paper lookup table convention. It is assumed 
that the detector to image correction is the same for all chips but it can be extended 
to arbitrary number of chips and extensions if necessary.

For ACS the correction is stored as an image extension with one row. Each element in 
the row specifies the correction in pixels for every pixel in the column (or row) in 
the science extension as predetermined by the calibration teams who would be responsible 
for creating the reference files. For ACS the correction is in the X direction and for 
WFPC2 - in the Y direction. The following new keywords are added to the header of each 
science extension of a science file: 

::

 'D2IMFILE' = "string - name of reference file to be used for creating the lookup table"
 'AXISCORR' = "integer (1 or 2) - axis to which the det2im correction is applied"
 'D2IMEXT' = "string - name of reference file which was last used to create the lookup table"
 'D2IMERR' = (optional)" float - maximum value of the correction"

'D2IMFILE' is used by UPDATEWCS as a flag that a reference file with this correction exists 
and an extension should be created. UPDATEWCS records the name of the reference file used 
for the lookup table extension to a keyword D2IMEXT in the primary header. It also populates 
keyword 'AXISCORR' based on whether this is a row or column correction. The lookup table 
extension has an 'EXTNAME' value of 'D2IMARR'.

'AXISCORR' is used as an indication of the axis to which the correction should be applied 
(1 - 'X' Axis, 2- 'Y' axis). 'D2IMEXT' stores the name of the reference file used by 
UPDATEWCS to create a ``D2IMARR`` extension. If 'D2IMEXT' is present in the 'SCI' extension 
header and is different from the current value of D2IMFILe in the primary header, the 
correction array in ``D2IMARR`` is updated. The optional keyword 'D2IMERR' allows a user to 
ignore this correction without modifying other header keywords by passing a parameter to 
the software. The HSTWCS class accepts a parameter 'minerr' which specifies the minimum 
value a distortion correction must have in order to be applied. If 'minerr' is larger than 
'D2IMERR' the correction is not applied. 

Detector To Image Reference File
--------------------------------

An entirely new reference file, the D2IMFILE reference table, serves as the source of this 1-D correction 
for each affected instrument. This reference file only contains a single array of offsets 
corresponding to the 1-D correction to be applied. Header keywords in the reference file 
specify which axis needs this correction. As a result, this new reference file remains 
small enough to easily be added to an input image without significant change in size. An 
initial **D2IMFILE** for ACS has been generated for testing with a sample header provided in 
:ref:`Appendix3`. 

.. _figure2:

.. figure:: /images/d2im_bar.png
   :width: 95 %
   :alt: ACS/WFC F475W D2IMFILE corrections
   :align: center
   
   This figure illustrates the corrections included in the first 246 columns of 
   the ACS/WFC F475W D2IMFILE.

The WCS for this correction describes the extension as a 1-D image, even though it gets 
applied to a 2-D image. This keeps it clear that the same correction gets applied to 
all rows(columns) without interpolation. The header specifies which axis this correction 
applies to through the use of the AXISCORR keyword. The WCS keywords in the header of the 
``D2IMARR`` extension specifies the transformation between pixel coordinates and lookup table 
position as if the lookup table were an image itself with 1-based positions (starting pixel 
is at a position of (1,1)). The value at that lookup table position then gets used to correct 
the original input pixel position.

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

Updating the FITS File
======================
Updating each science image with the distortion model using this merged
convention requires integrating these new reference files directly into the FITS file. 
This update gets performed using the following steps:

* determining what reference files should be applied to the science image
* read in distortion coefficients from IDCTAB reference file 
* [for ACS data only] compute time-dependent (TDD) skew terms from model described in IDCTAB file
* read in velocity aberration correction factor (VAFACTOR) keyword
* apply velocity aberration, and the TDD terms for ACS data as well, to the distortion coefficients

  * write time-corrected distortion coefficients as the SIP keywords

* [if d2imfile is to be applied] read in D2IMFILE reference table

  * update D2IMEXT with name of reference table and AXISCORR keyword with axis to be corrected
  * append D2IMFILE array as a new ``D2IMARR`` extension 

* [if NPOLFILE is to be applied] divide the NPOLFILE arrays by the linear distortion coefficients

  * write out normalized NPOLFILE arrays as new ``WCSDVARR`` extensions
  * update each SCI extension in the science image with the record-value keywords to point to the 2 ``WCSDVARR`` extensions (one for X corrections, one for Y corrections) associated with the SCI extension's chip

The STWCS task **updatewcs** applies these steps to update a science image's FITS file to 
incorporate the distortion model components using this convention. It not only modifies
the input reference file data to apply to each image to account for time-dependent and
velocity-aberration corrections as needed, but also creates the new extensions which get
appended to the science image's FITS file. 

Creating the D2IMARR extension
------------------------------
Converting the D2IMFILE reference table into a new ``D2IMARR`` FITS image extension involves only a few simple revisions 
to the header from D2IMFILE.  The header of the ``D2IMARR`` extension consists of the following keywords required in order to 
properly interpret and apply the data in the extension to the science array:

* AXISCORR : Direction in which the det2im correction is applied
* EXTNAME  : Set to 'D2IMARR'
* EXTVER   : Set to 1
* NAXIS    : Number of axes
* ``NAXISj`` : Size of each axis
* ``CRPIXj`` : Reference point for each axis, set at axis center
* ``CRVALj`` : computed from input science image array center on chip 
* ``CDELTj`` : Binning of axis, computed as :math:`1/BINAXIS_i` keyword from science image

These keywords supplement the standard FITS required keywords for an image extension, including such keywords as PCOUNT, GCOUNT, BITPIX, and XTENSION.
  
The corrections specified in this extension refer to pixel positions on the detector.  Since science images can be taken both as subarrays and in binned modes for some instruments, the subarray offset and binning factor get used to  compute the translation from science image pixel position into unbinned full-detector pixel positions.  Subarray exposures taken by HST detectors record the position of the detector's origin, (0,0) pixel, as ``LTVj`` keywords to identify what pixels on the physical detector were read out for the exposure. The conversion factor from image pixel position to physical detector pixel position of ``(NAXISj/2 + LTVj)*BINAXISj`` gets recorded as the ``CRVALj`` keyword value and gets used to correctly apply this correction to the science image. 

In addition to the pixel position transformations encoded as the ``D2IMARR`` WCS, keywords reporting how the D2IM correction was created get copied into the new ``D2IMARR`` image extension header from the primary header of the D2IMFILE.  This maintains as much provenance as possible for this correction. 

A full listing of the ``D2IMARR`` extension for a sample ACS image can be found in :ref:`d2imarr-header` in :ref:`Appendix1`. 


Creating the WCSDVARR Extension
-------------------------------
The NPOLFILE reference file contains at least 2 image extensions, one for the X correction and one for the Y correction for each chip. All these extensions get converted into their own ``WCSDVARR`` extension based on the FITS Distortion Paper convention when the NPOLFILE gets incorporated into the science image as another component of the distortion model. Both the array data for each NPOLFILE extension and the corresponding header needs to be modified before it can be written into the science image FITS file as a new ``WCSDVARR`` image extension. 

The data from the NPOLFILE arrays represent the residuals after accounting for the distortion model, yet this correction gets applied as part of the distortion correction described in :ref:`Equation 4 <equation4>`.  The linear terms of the distortion model need to be removed from the data in each NPOLFILE array in order to avoid applying the linear terms twice when applying the correction to the science data. This gets performed by reading in the linear distortion coefficients directly from the OCX and OCY keywords written out along with the SIP keywords, the multiplying them into the NPOLFILE data values using matrix dot operator to get the final, image specific NPOL correction to be written out as the ``WCSDVARR`` extension.

The header of this new ``WCSDVARR`` extension provides the translation from science image pixels to NPOLFILE array pixel positions as well as reporting on the provenance of the calibrations as recorded in the original NPOLFILE.  The following keywords get computed based on the values directly from the NPOLFILE header:

* ``NAXISj``  : Length of each axis
* ``CDELTj``  : Step size in detector pixels along each axis for the NPOL array
* ``CRPIXj``  : Reference pixel position of NPOL array
* ``CRVALj``  : Reference pixel position of NPOL array relative to science array
* EXTNAME          : always set to WCSDVARR
* EXTVER           : identifier reported in the DP.EXTVER record-value keywords in the science array header

These keywords supplement the standard FITS required keywords for an image extension, including such keywords as PCOUNT, GCOUNT, BITPIX, and XTENSION.  In addition, all keywords from the NPOLFILE primary header after and including 'FILENAME' get copied into the header of each WCSDARR extension to preserve the provenance of the calibration.  

The look-up tables are saved as separate FITS image extensions in the science files with ``EXTNAME`` 
set to ``WCSDVARR``. ``EXTVER`` is used when more than one look-up table is present in a single science 
file. Software which performs coordinate transformation will use bilinear interpolation to get 
the value of the distortion at a certain location in the image array. To fully map the image 
array to the distortion array the standard WCS keywords ``CRPIXj``, ``CRVALj`` and ``CDELTj`` are used. The 
mapping follows the transformation 

.. math:: 
   :label: Equation 6

    p_{j} = s_{j}(p_{j}-r_{j}) + w_{j}

where :math:`r_{j}` is the ``CRPIXj`` value in the distortion array which
corresponds to the :math:`w_{j}` value in the image array, recorded as
``CRVALj`` in the ``WCSDVARR`` header. Elements in the distortion array are spaced
by :math:`s_j` pixels in the image array, where :math:`s_j` is the ``CDELTj``
value in the distortion array header.  In general :math:`s_j` can have
a non-integer value but cannot be zero. However, if the distortion array
was obtained as a subimage of a larger array having a non-integer step size
can produce undesirable results during interpolation. A full listing of the 
``WCSDVARR`` extension for a sample ACS image can be found in :ref:`wcsdvarr-header` in :ref:`Appendix1`. 

Summary
=======
This paper describes a merging of previously proposed FITS WCS conventions to fully support the multi-component distortion models derived from calibrations for HST detectors.  The application of this merged convention allows each science image to contain the full distortion model applicable to that specific image in an efficient and FITS compatible manner.  The use of this calibration in the DrizzlePac package has been demonstrated to correct science data to much better than 0.1 pixels across each image's field of view, with a typical RMS for aligning two ACS images on the order of 0.03 pixels in a suitably dense field of sources. This convention, despite making a few basic assumptions, retains each separate FITS convention's full functionality so that any software which understood, for example, the SIP standard will still work as before with the SIP keywords written out by the convention.  

All HST ACS and WFC3 images retrieved from the archive have been updated using this convention so that users will no longer need to retrieve the distortion calibration data separately. Anyone using HST images will now be able to use the STWCS and/or DrizzlePac package to perform coordinate transformations or image alignment based on this convention, while still being able to use external tools like DS9 to take advantage of the SIP conventions as well. This solution now provides the best possible solution for supporting these highly accurate, yet complex multi-component distortion models in the most efficient manner available to data written out in the FITS format. 

References
==========
.. [DistortionPaper] Calabretta M. R., Valdes F. G., Greisen E. W., and Allen S. L., 2004, 
    "Representations of distortions in FITS world coordinate systems",[cited 2012 Sept 18], 
    Available from: http://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf

.. [SIPConvention] Shupe D.L., Hook R.N., 2008, "The SIP Convention for Representing Distortion in FITS Image
    Headers", [cited 2012 Sept 18], Available from: http://fits.gsfc.nasa.gov/registry/sip.html


.. _Appendix1:

==================================
Appendix 1 - Sample ACS/WFC Image 
==================================

The WCS of a single chip from an ACS/WFC exposure illustrates how the SIP keywords are
populated based on the coefficients from the external IDCTAB reference file.  In addition,
this header includes the keywords referring to additional distortion corrections
related to non-polynomial corrections from the NPOLFILE and to column-width corrections from
the D2IMFILE.  This sample illustrates how all three corrections can be specified at the
same time in a FITS header using our rules for combining the SIP WCS convention and
FITS Distortion Paper proposed syntax, while also using FITS WCS Paper I alternate WCS 
standards to maintain a record of the WCS information prior to being updated/recomputed to
use the new reference information. The old WCS gets stored using WCS key 'O' and 'WCSNAMEO' = 'OPUS'
to indicate it was originally computed by OPUS, the HST pipeline system. 

FITS File extensions
==================================

The FITS file for this ACS/WFC image now contains extra extensions for the NPOLFILE and D2IMFILE
corrections.

::

 Filename: jbf401p8q_flc.fits
 No.    Name         Type      Cards   Dimensions   Format
 0    PRIMARY     PrimaryHDU     261   ()           int16   
 1    SCI         ImageHDU       184   (4096, 2048)   float32   
 2    ERR         ImageHDU        55   (4096, 2048)   float32   
 3    DQ          ImageHDU        47   (4096, 2048)   int16   
 4    SCI         ImageHDU       183   (4096, 2048)   float32   
 5    ERR         ImageHDU        55   (4096, 2048)   float32   
 6    DQ          ImageHDU        47   (4096, 2048)   int16   
 7    D2IMARR     ImageHDU        12   (4096,)      float32   
 8    WCSDVARR    ImageHDU        37   (65, 33)     float32   
 9    WCSDVARR    ImageHDU        37   (65, 33)     float32   
 10   WCSDVARR    ImageHDU        37   (65, 33)     float32   
 11   WCSDVARR    ImageHDU        37   (65, 33)     float32   
 12   WCSCORR     BinTableHDU     59   14R x 24C    [40A, I, 1A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 
                                                     24A, 24A, D, D, D, D, J, 40A, 128A] 

The last extension, named WCSCORR, contains a binary table providing a summary of all the WCS 
solutions that have been applied to this file and does not act as an active part of the WCS 
or its interpretation.

Primary Header
==================================

The PRIMARY header of HST data contains keywords specifying information general to 
the entire file, such as what calibration steps were applied and what reference files
should be used.  No active WCS keywords (keywords interpreted for coordinate transformations)
are present in the PRIMARY header, but keywords specifying the applicable distortion 
reference files can be found in the PRIMARY header. Some keywords describing the 
distortion model and when the WCS was updated with the distortion information from the
reference files can also be found in the PRIMARY header. These distortion and WCS
related keywords from the PRIMARY header are::


              / CALIBRATION REFERENCE FILES                                     
                                                                                
 IDCTAB  = 'jref$v8q1444sj_idc.fits' / image distortion correction table         
 DGEOFILE= 'jref$qbu16420j_dxy.fits' / Distortion correction image               
 D2IMFILE= 'jref$v971826mj_d2i.fits' / Column Correction Reference File          
 NPOLFILE= 'jref$v971826aj_npl.fits' / Non-polynomial Offsets Reference File     

 UPWCSVER= '1.0.0   '           / Version of STWCS used to updated the WCS       
 PYWCSVER= '1.11-4.10'          / Version of PYWCS used to updated the WCS       
 DISTNAME= 'jbf401p8q_v8q1444sj-v971826aj-v971826mj'                             
 SIPNAME = 'jbf401p8q_v8q1444sj'                                                 

The remainder of the PRIMARY header specifies the general characteristics of the image.


SCI Header Keywords
==================================

The following keywords only represent the WCS keywords from a sample ACS/WFC SCI header with 4-th order
polynomial distortion correction from the IDCTAB reference file, along with NPOLFILE and 
D2IMFILE corrections from the specific reference files used as examples in :ref:`Appendix2`
:ref:`Appendix3`.

::

 XTENSION= 'IMAGE   '           / IMAGE extension                                
 BITPIX  =                  -32                                                  
 NAXIS   =                    2                                                  
 NAXIS1  =                 4096                                                  
 NAXIS2  =                 2048                                                  
 PCOUNT  =                    0 / required keyword; must = 0                     
 GCOUNT  =                    1 / required keyword; must = 1                     
 ORIGIN  = 'HSTIO/CFITSIO March 2010'                                            
 DATE    = '2012-06-13' / date this file was written (yyyy-mm-dd)                
 INHERIT =                    T / inherit the primary header                     
 EXTNAME = 'SCI     '           / extension name                                 
 EXTVER  =                    1 / extension version number                       
 ROOTNAME= 'jbf401p8q                         ' / rootname of the observation set
 EXPNAME = 'jbf401p8q                ' / exposure identifier                     
 BUNIT   = 'ELECTRONS'          / brightness units                               
                                                                                
              / WFC CCD CHIP IDENTIFICATION                                     
                                                                                
 CCDCHIP =                    2 / CCD chip (1 or 2)                              
                                                                                
              / World Coordinate System and Related Parameters                  
                                                                                
 WCSAXES =                    2 / number of World Coordinate System axes         
 CRPIX1  =                 2048 / x-coordinate of reference pixel                
 CRPIX2  =                 1024 / y-coordinate of reference pixel                
 CRVAL1  =        11.3139376926 / first axis value at reference pixel            
 CRVAL2  =        42.0159325283 / second axis value at reference pixel           
 CTYPE1  = 'RA---TAN-SIP'       / the coordinate type for the first axis         
 CTYPE2  = 'DEC--TAN-SIP'       / the coordinate type for the second axis        
 CD1_1   = -7.8194868997837E-06 / partial of first axis coordinate w.r.t. x      
 CD1_2   = 1.09620231564470E-05 / partial of first axis coordinate w.r.t. y      
 CD2_1   = 1.14279318521882E-05 / partial of second axis coordinate w.r.t. x     
 CD2_2   = 8.66885775536641E-06 / partial of second axis coordinate w.r.t. y     
 LTV1    =        0.0000000E+00 / offset in X to subsection start                
 LTV2    =        0.0000000E+00 / offset in Y to subsection start                
 LTM1_1  =                  1.0 / reciprocal of sampling rate in X               
 LTM2_2  =                  1.0 / reciprocal of sampling rate in Y               
 ORIENTAT=    51.66276166150634 / position angle of image y axis (deg. e of n)   
 RA_APER =   1.133205840898E+01 / RA of aperture reference position              
 DEC_APER=   4.202747924810E+01 / Declination of aperture reference position     
 PA_APER =              51.4653 / Position Angle of reference aperture center (de
 VAFACTOR=   9.999374411935E-01 / velocity aberration plate scale factor         
                                                                                 
 WCSCDATE= '18:41:12 (13/06/2012)' / Time WCS keywords were copied.              
 A_0_2   = 2.18045745103211E-06                                                  
 B_0_2   = -7.2266555836441E-06                                                  
 A_1_1   = -5.2225148886672E-06                                                  
 B_1_1   = 6.20296011911662E-06                                                  
 A_2_0   = 8.54842918202735E-06                                                  
 B_2_0   = -1.7551668097547E-06                                                  
 A_0_3   = 8.09354090167772E-12                                                  
 B_0_3   = -4.2488740853874E-10                                                  
 A_1_2   = -5.2903025382457E-10                                                  
 B_1_2   = -7.6098727022982E-11                                                  
 A_2_1   = -4.4821374838034E-11                                                  
 B_2_1   = -5.1244088812978E-10                                                  
 A_3_0   = -4.6755353102513E-10                                                  
 B_3_0   = 8.48145748580355E-11                                                  
 A_0_4   = -8.3665541956904E-17                                                  
 B_0_4   = -2.1662072760964E-14                                                  
 A_1_3   = -1.5108585176304E-14                                                  
 B_1_3   = -1.5686763638364E-14                                                  
 A_2_2   = 3.61252682019403E-14                                                  
 B_2_2   = -2.6194214315839E-14                                                  
 A_3_1   = 1.03502537140899E-14                                                  
 B_3_1   = -2.6915637616404E-15                                                  
 A_4_0   = 2.32643027828425E-14                                                  
 B_4_0   = -1.5701287138447E-14                                                  
 A_ORDER =                    4                                                  
 B_ORDER =                    4                                                  
 IDCSCALE=                 0.05                                                  
 IDCV2REF=    256.6019897460938                                                  
 IDCV3REF=    302.2520141601562                                                  
 IDCTHETA=                  0.0                                                  
 OCX10   = 0.001965125839177266                                                  
 OCX11   =  0.04983026381230307                                                  
 OCY10   =   0.0502766128737329                                                  
 OCY11   = 0.001493971240339153                                                  
 TDDALPHA=    0.246034678162242                                                  
 TDDBETA = -0.07934489272074734                                                  
 IDCXREF =               2048.0                                                  
 IDCYREF =               1024.0                                                  
 AXISCORR=                    1                                                  
 D2IMEXT = '/grp/hst/cdbs/jref/v971826mj_d2i.fits'                               
 D2IMERR = 0.002770500956103206                                                  
 WCSNAMEO= 'OPUS    '                                                            
 WCSAXESO=                    2                                                  
 CRPIX1O =                 2048                                                  
 CRPIX2O =                 1024                                                  
 CDELT1O =                    1                                                  
 CDELT2O =                    1                                                  
 CUNIT1O = 'deg     '                                                            
 CUNIT2O = 'deg     '                                                            
 CTYPE1O = 'RA---TAN-SIP'                                                        
 CTYPE2O = 'DEC--TAN-SIP'                                                        
 CRVAL1O =        11.3139376926                                                  
 CRVAL2O =        42.0159325283                                                  
 LONPOLEO=                  180                                                  
 LATPOLEO=        42.0159325283                                                  
 RESTFRQO=                    0                                                  
 RESTWAVO=                    0                                                  
 CD1_1O  =   -7.81948731152E-06                                                  
 CD1_2O  =    1.09620228331E-05                                                  
 CD2_1O  =    1.14279315609E-05                                                  
 CD2_2O  =    8.66885813904E-06                                                  
 WCSNAME = 'IDC_v8q1444sj'                                                       
 CPERR1  =                  0.0 / Maximum error of NPOL correction for axis 1    
 CPDIS1  = 'Lookup  '           / Prior distortion funcion type                  
 DP1     = 'EXTVER: 1' / Version number of WCSDVARR extension containing lookup d
 DP1     = 'NAXES: 2' / Number of independent variables in distortion function   
 DP1     = 'AXIS.1: 1' / Axis number of the jth independent variable in a distort
 DP1     = 'AXIS.2: 2' / Axis number of the jth independent variable in a distort
 CPERR2  =                  0.0 / Maximum error of NPOL correction for axis 2    
 CPDIS2  = 'Lookup  '           / Prior distortion funcion type                  
 DP2     = 'EXTVER: 2' / Version number of WCSDVARR extension containing lookup d
 DP2     = 'NAXES: 2' / Number of independent variables in distortion function   
 DP2     = 'AXIS.1: 1' / Axis number of the jth independent variable in a distort
 DP2     = 'AXIS.2: 2' / Axis number of the jth independent variable in a distort
 NPOLEXT = 'jref$v971826aj_npl.fits'                                             


All keywords related to the exposure itself, such as readout pattern, have been deleted 
from this SCI header listing for the sake of brevity. 

.. _d2imarr-header:

D2IMARR Header
==============
The full, complete header of the ``D2IMARR`` extension as derived from the D2IMFILE 
discussed in :ref:`Appendix3`.

::

 XTENSION= 'IMAGE   '           / Image extension                                
 BITPIX  =                  -32 / array data type                                
 NAXIS   =                    1 / number of array dimensions                     
 NAXIS1  =                 4096                                                  
 PCOUNT  =                    0 / number of parameters                           
 GCOUNT  =                    1 / number of groups                               
 AXISCORR=                    1 / Direction in which the det2im correction is app
 EXTVER  =                    1 / Distortion array version number                
 EXTNAME = 'D2IMARR '           / WCS distortion array                           
 CDELT1  =                  1.0 / Coordinate increment along axis                
 CRPIX1  =               2048.0 / Coordinate system reference pixel              
 CRVAL1  =               2048.0 / Coordinate system value at reference pixel     

.. _wcsdvarr-header:

WCSDVARR Header
===============
Each of the WCSDVARR extensions has been derived based on the values for the 
NPOL correction found in the reference file described in :ref:`Appendix2`. The 
full header for the WCSDVARR extension with EXTVER=1 is::

 XTENSION= 'IMAGE   '           / Image extension                                
 BITPIX  =                  -32 / array data type                                
 NAXIS   =                    2 / number of array dimensions                     
 NAXIS1  =                   65                                                  
 NAXIS2  =                   33                                                  
 PCOUNT  =                    0 / number of parameters                           
 GCOUNT  =                    1 / number of groups                               
 EXTVER  =                    1 / Distortion array version number                
 EXTNAME = 'WCSDVARR'           / WCS distortion array                           
 CRVAL2  =                  0.0 / Coordinate system value at reference pixel     
 CRPIX1  =                  0.0 / Coordinate system reference pixel              
 CRPIX2  =                  0.0 / Coordinate system reference pixel              
 CRVAL1  =                  0.0 / Coordinate system value at reference pixel     
 CDELT1  =                   64 / Coordinate increment along axis                
 CDELT2  =                   64 / Coordinate increment along axis                
 FILENAME= 'v971826aj_npl.fits' / name of file                                   
 FILETYPE= 'DXY GRID'           / type of data found in data file                
 OBSTYPE = 'IMAGING '           / type of observation                            
 TELESCOP= 'HST'                / telescope used to acquire data                 
 INSTRUME= 'ACS   '             / identifier for instrument used to acquire data 
 DETECTOR= 'WFC'                / detector in use: WFC, HRC, or SBC              
 FILTER1 = 'F475W   '           / element selected from filter wheel 1           
 FILTER2 = 'CLEAR2L '           / element selected from filter wheel 2           
 USEAFTER= 'Mar 01 2002 00:00:00'                                                
 COMMENT = 'NPOL calibration file created by Ray A. Lucas 29 APR 2010'           
 DESCRIP = 'Residual geometric distortion file for use with astrodrizzle-------' 
 PEDIGREE= 'INFLIGHT 11/11/2002 11/11/2002'                                      
 HISTORY   Non-polynomial offset file generated from qbu16420j_dxy.fits          
 HISTORY   Only added to the flt.fits file and used in coordinate                
 HISTORY   transformations if the npol reference filename is specified in        
 HISTORY   the header.  The offsets are copied from the reference file into      
 HISTORY   two arrays for each chip.  Each array is stored as a 65x33 pixel      
 HISTORY   image that gets interpolated up to the full chip size. Two new        
 HISTORY   extensions for each chip are also appended to the flt file            
 HISTORY   (WCSDVARR).                                                           
 HISTORY qbu16420j_npl.fits renamed to v9615069j_npl.fits on Sep 6 2011          
 HISTORY v9615069j_npl.fits renamed to v971826aj_npl.fits on Sep 7 2011 
 
Each of the ``WCSDVARR`` extension headers contains the same set of keywords, with
only the values varying to reflect the axis and chip corrected by this extension.




.. _Appendix2:

==================================
Appendix 2 - NPOLFILE Example 
==================================
The NPOLFILE reference file format includes a PRIMARY header describing what kind of 
image should be corrected by this file, along with extensions containing the corrections
for each chip.  

FITS File Extensions
==================================

A sample NPOLFILE applicable to ACS/WFC F475W images has the FITS extensions::

 Filename: /grp/hst/cdbs/jref/v971826aj_npl.fits
 No.    Name         Type      Cards   Dimensions   Format
 0    PRIMARY     PrimaryHDU      35   ()           int16   
 1    DX          ImageHDU       180   (65, 33)     float32   
 2    DY          ImageHDU       215   (65, 33)     float32   
 3    DX          ImageHDU       215   (65, 33)     float32   
 4    DY          ImageHDU       215   (65, 33)     float32   

The extensions with the name 'DX' provide the corrections in X for each of the 
ACS/WFC's 2 chips, while the 'DY' extensions provide the corrections in Y for each chip.

Primary Header
==================================

The PRIMARY header of this file only includes the minimum information necessary to describe
what exposures should be corrected by this reference file and how it was generated. A full
listing of the PRIMARY header includes::

 SIMPLE  =                    T / Fits standard                                  
 BITPIX  =                   16 / Bits per pixel                                 
 NAXIS   =                    0 / Number of axes                                 
 EXTEND  =                    T / File may contain extensions                    
 ORIGIN  = 'NOAO-IRAF FITS Image Kernel July 2003' / FITS file originator        
 IRAF-TLM= '2011-09-09T13:24:40'                                                 
 NEXTEND =                    4 / Number of standard extensions                  
 DATE    = '2010-04-02T19:53:08'                                                 
 FILENAME= 'v971826aj_npl.fits' / name of file                                   
 FILETYPE= 'DXY GRID'           / type of data found in data file                
 OBSTYPE = 'IMAGING '           / type of observation                            
 TELESCOP= 'HST'                / telescope used to acquire data                 
 INSTRUME= 'ACS   '             / identifier for instrument used to acquire data 
 DETECTOR= 'WFC'                / detector in use: WFC, HRC, or SBC              
 FILTER1 = 'F475W   '           / element selected from filter wheel 1           
 FILTER2 = 'CLEAR2L '           / element selected from filter wheel 2           
 USEAFTER= 'Mar 01 2002 00:00:00'                                                
 COMMENT = 'NPOL calibration file created by Ray A. Lucas 29 APR 2010'           
 DESCRIP = 'Residual geometric distortion file for use with astrodrizzle-------' 
 PEDIGREE= 'INFLIGHT 11/11/2002 11/11/2002'                                      
 HISTORY   Non-polynomial offset file generated from qbu16420j_dxy.fits          
 HISTORY   Only added to the flt.fits file and used in coordinate                
 HISTORY   transformations if the npol reference filename is specified in        
 HISTORY   the header.  The offsets are copied from the reference file into      
 HISTORY   two arrays for each chip.  Each array is stored as a 65x33 pixel      
 HISTORY   image that gets interpolated up to the full chip size. Two new        
 HISTORY   extensions for each chip are also appended to the flt file            
 HISTORY   (WCSDVARR).                                                           
 HISTORY qbu16420j_npl.fits renamed to v9615069j_npl.fits on Sep 6 2011          
 HISTORY v9615069j_npl.fits renamed to v971826aj_npl.fits on Sep 7 2011          


Data Extension Header
==================================

Each ACS/WFC chip has a shape of 4096 x 2048 pixels,
yet the data arrays in this specific reference file only have 65x33 values.
Each data extension ('DX' and 'DY') contains only those keywords necessary to 
properly interpolate the sub-sampled values from the arrays to apply to each individual
pixel in the full ACS/WFC exposure. The full header for the ['DX',1] extension contains::

 XTENSION= 'IMAGE   '           / Image extension                                
 BITPIX  =                  -32 / Bits per pixel                                 
 NAXIS   =                    2 / Number of axes                                 
 NAXIS1  =                   65 / Axis length                                    
 NAXIS2  =                   33 / Axis length                                    
 PCOUNT  =                    0 / No 'random' parameters                         
 GCOUNT  =                    1 / Only one group                                 
 EXTNAME = 'DX      '           / Extension name                                 
 EXTVER  =                    1 / Extension version                              
 ORIGIN  = 'NOAO-IRAF FITS Image Kernel July 2003' / FITS file originator        
 INHERIT =                    F / Inherits global header                         
 DATE    = '2004-04-28T16:44:21'                                                 
 IRAF-TLM= '16:42:00 (30/11/2006)'                                               
 WCSDIM  =                    2                                                  
 LTM1_1  =                   1.                                                  
 LTM2_2  =                   1.                                                  
 WAT0_001= 'system=physical'                                                     
 WAT1_001= 'wtype=linear'                                                        
 WAT2_001= 'wtype=linear'                                                        
 CCDCHIP =                    2 / CCDCHIP from full size dgeo file               
 LTV1    =                    0                                                  
 LTV2    =                    0                                                  
 ONAXIS1 =                 4096 / NAXIS1 of full size dgeo file                  
 ONAXIS2 =                 2048 / NAXIS2 of full size dgeo file                  
 CDELT1  =                   64 / Coordinate increment along axis                
 CDELT2  =                   64 / Coordinate increment along axis                


.. _Appendix3:

==================================
Appendix 3 - D2IMFILE Example 
==================================

The D2IMFILE reference file only contains a single 1-D array that should correct the
column (row) values based on the value of the 'AXISCORR' keyword in the SCI header. 

FITS File Extensions
==================================

This simple reference file, therefore, contains only 2 extensions; namely,

::

 Filename: /grp/hst/cdbs/jref/v971826mj_d2i.fits
 No.    Name         Type      Cards   Dimensions   Format
 0    PRIMARY     PrimaryHDU      35   ()           int16   
 1    DX          ImageHDU        18   (4096,)      float32   

PRIMARY Header
==================================

The PRIMARY header only needs to contain information on what detector this file corrects,
along with any available information on how this file was generated.  The ACS/WFC D2IMFILE
PRIMARY header only includes::

 SIMPLE  =                    T / Fits standard                                  
 BITPIX  =                   16 / Bits per pixel                                 
 NAXIS   =                    0 / Number of axes                                 
 EXTEND  =                    T / File may contain extensions                    
 ORIGIN  = 'NOAO-IRAF FITS Image Kernel July 2003' / FITS file originator        
 DATE    = '2010-02-01T20:19:11' / Date FITS file was generated                  
 IRAF-TLM= '2011-09-02T13:04:07' / Time of last modification                     
 NEXTEND =                    1 / number of extensions in file                   
 FILENAME= 'v971826mj_d2i.fits' / name of file                                   
 FILETYPE= 'WFC D2I FILE'          / type of data found in data file             
 OBSTYPE = 'IMAGING '              / type of observation                         
 TELESCOP= 'HST'                / telescope used to acquire data                 
 INSTRUME= 'ACS   '             / identifier for instrument used to acquire data 
 DETECTOR= 'WFC     '                                                            
 USEAFTER= 'Mar 01 2002 00:00:00'                                                
 COMMENT = 'D2I calibration file created by Warren Hack 29 APR 2010'             
 DESCRIP = 'Column-width correction file for WFC images------------------------' 
 PEDIGREE= 'INFLIGHT 11/11/2002 11/11/2002'                                      
 HISTORY                                                                         
 HISTORY   Fixed column (or row) width correction file. This is applied          
 HISTORY   as a correction to the input pixel position and the output of         
 HISTORY   this correction is to be used as input to the polynomial and          
 HISTORY   non-polynomial distortion corrections.                                
 HISTORY                                                                         
 HISTORY   For ACS WFC data, the correction is stored as an image extension      
 HISTORY   (D2IMARR) with one row. Each element in the row specifies the         
 HISTORY   correction in pixels for every pixel in the column (or row) in        
 HISTORY   the science extension; for ACS WFC, the correction is in the X        
 HISTORY   direction.                                                            
 HISTORY                                                                         
 HISTORY   For a more in-depth explanation of this file, please see the          
 HISTORY   draft writeup at:                                                     
 HISTORY http://stsdas.stsci.edu/stsci_python_epydoc/stwcs/fits_conventions.html 
 HISTORY wfc_ref68col_d2i.fits renamed to v961506lj_d2i.fits on Sep 6 2011       
 HISTORY v961506lj_d2i.fits renamed to v971826mj_d2i.fits on Sep 7 2011          

In this case, most of the keywords not required by FITS describe how this file
was computed while also describing how it should be applied. 

Data Extension Header
==================================

The header keywords for the actual DX array simply needs to provide the information
necessary to apply the values to the data; namely, 

::

 XTENSION= 'IMAGE   '           / Image extension                                
 BITPIX  =                  -32 / Bits per pixel                                 
 NAXIS   =                    1 / Number of axes                                 
 NAXIS1  =                 4096 / Axis length                                    
 PCOUNT  =                    0 / No 'random' parameters                         
 GCOUNT  =                    1 / Only one group                                 
 EXTNAME = 'DX      '           / Extension name                                 
 EXTVER  =                   11 / Extension version                              
 ORIGIN  = 'NOAO-IRAF FITS Image Kernel July 2003' / FITS file originator        
 INHERIT =                    F / Inherits global header                         
 DATE    = '2009-03-18T19:28:09' / Date FITS file was generated                  
 IRAF-TLM= '16:05:02 (18/03/2009)' / Time of last modification                   
 CRPIX1  =                    0 / Distortion array reference pixel               
 CDELT1  =                    0 / Grid step size in first coordinate             
 CRVAL1  =                    0 / Image array pixel coordinate                   
 CRPIX2  =                    0 / Distortion array reference pixel               
 CDELT2  =                    0 / Grid step size in second coordinate            
 CRVAL2  =                    0 / Image array pixel coordinate                   

The fact that these values get applied without interpolation to each pixel in a row,
in this case, means that no translation terms are needed in the header, making for 
a very simple header and very simple application to the data.
