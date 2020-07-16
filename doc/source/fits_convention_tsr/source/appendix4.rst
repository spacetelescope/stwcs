.. _appendix4:

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
--------------------

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
---------------

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
--------------------

The following keywords only represent the WCS keywords from a sample ACS/WFC SCI header with 4-th order
polynomial distortion correction from the IDCTAB reference file, along with NPOLFILE and
D2IMFILE corrections from the specific reference files used as examples in :ref:`appendix2`
:ref:`appendix3`.

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
 DP1     = 'EXTVER: 1' / Version number of WCSDVARR extension
 DP1     = 'NAXES: 2' / Number of independent variables in CPDIS function
 DP1     = 'AXIS.1: 1' / Axis number of the 1st variable in a CPDIS function
 DP1     = 'AXIS.2: 2' / Axis number of the 2nd variable in a CPDIS function
 CPERR2  =                  0.0 / Maximum error of NPOL correction for axis 2
 CPDIS2  = 'Lookup  '           / Prior distortion funcion type
 DP2     = 'EXTVER: 2' / Version number of WCSDVARR extension
 DP2     = 'NAXES: 2' / Number of independent variables in CPDIS function
 DP2     = 'AXIS.1: 1' / Axis number of the 1st variable in a CPDIS function
 DP2     = 'AXIS.2: 2' / Axis number of the 2nd variable in a CPDIS function
 NPOLEXT = 'jref$v971826aj_npl.fits'


All keywords related to the exposure itself, such as readout pattern, have been deleted
from this SCI header listing for the sake of brevity.

.. _d2imarr-header:

D2IMARR Header
--------------------

The full, complete header of the ``D2IMARR`` extension as derived from the D2IMFILE
discussed in :ref:`appendix3`.

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
--------------------

Each of the WCSDVARR extensions has been derived based on the values for the
NPOL correction found in the reference file described in :ref:`appendix2`. The
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
