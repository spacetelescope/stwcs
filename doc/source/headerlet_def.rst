===========================================================================
Definition of a Headerlet and It's Role in Updating WCS Information - DRAFT
===========================================================================

.. abstract::
   :author: Warren Hack, Nadezhda Dencheva
   :date: 12 Oct 2010

   The 'headerlet' serves as a mechanism for encapsulating WCS information
   which can be used to update the WCS solution of an image. This object
   needs to be as compact as possible while providing an unambigious and
   self-consistent WCS solution for an image while requiring a minimum
   level of software necessary to apply the headerlet to an image.
   These basic requirements come from the desire to create a mechanism
   for passing improved astrometric solutions for HST data and provide
   those solutions in a manner that would not require getting entirely
   new images from an archive when only the WCS information has been
   updated. This report describes the format and contents of a headerlet
   along with the procedures which would be used to update images with
   the updated WCS information from the headerlet.

Introduction
============
The HST Archive provides access to all HST data, most of which can not be readily combined together due to errors in guide star astrometry imposing offsets between images taken using different pairs of guide stars.  A lot of effort has gone into computing those offsets so that all the images taken at a particular pointing can be combined successfully with astrometry matched to external astrometric catalogs. Unfortunately, there is no current mechanism for passing those updated solutions along to the community without providing entirely new copies of all the data.  

The concept of a 'headerlet' seeks to provide a solution where only the WCS solution for an image that has been aligned to an astrometric catalog can be archived and retrieved for use in updating copies of that image's WCS information without getting the image data again.  Multiple 'headerlets' could even be provided with each representing the alignment of an image to a different astrometric solution, giving the end user the option to get the solution that would allow them to best align their images with external data of interest to them.  These benefits can only be realized with the proper definition of a 'headerlet' and the procedures used to define them and apply them to data. 

Source Image
============
All users get their copy of an image from the HST Archive (OTFR) after being processed with the latest image calibrations, including applying the latest available distortion models. The calibrated image gets passed to the user with the ``_flt.fits`` suffix and is referred to as the ``FLT`` image.  These FLT images serve as the inputs to MultiDrizzle in order to apply the distortion models and combine the images into a single ``DRZ`` product.  

The WCS information in the FLT images has been updated to include the full distortion model, including the full polynomial solution from the IDCTAB and all the corrections formerly combined into the DGEOFILE. The :ref:`FITS Conventions Report <fits_conventions_tsr>` report by Dencheva contains the full description of the conventions used to describe all these components in a FITS file. The header now contains the following set of keywords and extensions to fully describe the WCS with distortion:

  * **Linear WCS keywords**: specifically, CRPIX, CRVAL, CTYPE, CD matrix keywords
  * **SIP coefficients**: A_*_* and B_*_*, A_ORDER, B_ORDER, 
    OCX10, OCX11, OCY10, and OCY11 keywords
  * **NPOL file**: if an NPOLFILE has been specified for the image, 
    CPDIS and DP record-value keywords to point to WCSDVARR extensions (FITS Distortion Paper convention)
  * **Column correction file**: if a D2IMFILE has been specified for use with the image, 
    the D2IMEXT and D2IMERR keywords signify the use of a D2IM file extension
  * **WCSDVARR extensions**: 2 extensions for each chip with lookup tables containing 
    the NPOLFILE corrections, with each extension corresponding to an axis of 
    the image (X correction or Y correction) (if an NPOLFILE has been applied to the image)
  * **D2IMVARR extension**: an extension with a lookup table containing the 
    column-correction for that chip from the D2IMFILE.
 

Each science header will have its own set of these keywords and extensions that will be kept together as part of the headerlet definition.  This avoids any ambiguity as to what solution was used for any given WCS. 

An HST ACS/WFC exposure would end up with the following set of extensions::

    EXT#  FITSNAME      FILENAME              EXTVE DIMENS       BITPI OBJECT       

    0     j8hw27c4q_flt j8hw27c4q_flt.fits                       16                 
    1       IMAGE       SCI                   1     4096x2048    -32                
    2       IMAGE       ERR                   1     4096x2048    -32                
    3       IMAGE       DQ                    1     4096x2048    16                 
    4       IMAGE       SCI                   2     4096x2048    -32                
    5       IMAGE       ERR                   2     4096x2048    -32                
    6       IMAGE       DQ                    2     4096x2048    16                 
    7       IMAGE       D2IMARR               1     4096         -32                
    8       IMAGE       WCSDVARR              1     64x32        -32                
    9       IMAGE       WCSDVARR              2     64x32        -32                
    10      IMAGE       WCSDVARR              3     64x32        -32                
    11      IMAGE       WCSDVARR              4     64x32        -32                

There may be a lot of extensions appended to this FITS file, but the sum total of all these new extensions comes to approximately 100kB for ACS/WFC images (our sample only requires 86400 bytes), making them a space efficient means of managing all the distortion and WCS information. 

Headerlet Definition
====================
The `headerlet` needs to be a self-consistent, fully described definition of a WCS and its distortion for all chips/detectors of a single exposure.  The WCS and SIP coefficients get derived from the SCI header directly, along with the all keywords which refer to the extensions with the optional NPOLFILE and D2IMFILE corrections.  The full 'headerlet' will be stored as a multi-extension FITS object that includes a primary header, an extension for each chip which only contains the WCS and distortion keywords, and any additional extensions for optional distortion correction information.  This object can be written out to a file and/or attached to an existing image's FITS file as a new extension.

The science observation's original WCS will be saved to a headerlet so that the user can always revert back to the original WCS solution at any time. 
There will be an option to permanently delete the original WCS and not save it to a headerlet.  

New WCS Extension
-----------------
A new extension, named SIPWCS, containing all the WCS-related keywords from the SCI header, including all keywords referring to NPOL and D2IM extensions as well as all sets of alternate WCS keywords, will be created to serve as the record of the original WCS. Keywords (TBD) recording the alignment information are recorded in this header as well. All the sets of linear WCS keywords stored using FITS Paper I Multiple WCS Standard would be defined using the same set of distortion coefficients written to the SIP keywords and NPOL files.  This insures that all the information in the header remains consistent. The keywords in this extension can be used to overwrite the keywords in the corresponding SCI header to update the WCS solution for each chip without any further modification or computation. The new extension then serves not only as a record of all the WCS solutions derived for the image, but also the source of values for restoring the SCI header WCS when desired.  


Headerlet File Structure
-----------------------------
This new extension along with the NPOLFILE and the D2IMFILE extensions fully describe the WCS of each chip and can serve without further modification as the definition of the `headerlet`. The listing of the FITS extensions for a `headerlet` for the sample ACS/WFC exposure after writing it out to a file would then be::

    EXT#  FITSNAME      FILENAME              EXTVE DIMENS       BITPI OBJECT       

    0     j8hw27c4q     j8hw27c4q_hdr.fits                       16
    1       IMAGE       SIPWCS                1                  8
    2       IMAGE       SIPWCS                2                  8
    3       IMAGE       WCSDVARR              1     64x32        -32                
    4       IMAGE       WCSDVARR              2     64x32        -32                
    5       IMAGE       WCSDVARR              3     64x32        -32                
    6       IMAGE       WCSDVARR              4     64x32        -32                
    7       IMAGE       D2IMARR               1     4096         -32                

This file now fully describes the WCS solution for this image, complete with all the distortion information used to originally define the solution. No further reference files or computations would be needed when this `headerlet` gets used to update an image.

.. note::

   A headerlet derived from a full-frame WFC3/UVIS image would only
   contain a PRIMARY header and two SIPWCS  extensions (one for each SCI extension)
   as WFC3/UVIS does not currently use NPOLFILE or D2IMFILE reference files as
   part of their distortion model.

The keywords used to populate the headerlet come from all the extensions of the updated
FITS file, as illustrated in the following figure.

.. figure:: images/Headerlet_figure_final.png
   :width: 95 %
   :alt: Relationship Between an ACS/WFC Image’s FITS File and a Headerlet
   :align: center

   This figure shows the keywords that are included in a headerlet, the extensions included 
   in a headerlet, and how a headerlet appears as a new extension when it gets appended to the
   original ACS/WFC file.


Headerlet Primary Header
-------------------------
The primary header must have 4 required keywords:

 * `HDRNAME`  - a unique name for the headerlet
 * `DESTIM`   - target image filename (the ROOTNAME keyword of the original archive filename)
 * `STWCSVER` - version of STWCS used to create the WCS of the original image
 * `PYWCSVER` - version of PyWCS used to create the WCS of the original image

These keywords are used for determining whether a headerlet can be applied to a
given exposure and how it needs to be applied. Additional keywords provide more
information about the solution itself, how it was derived, and by whom, through use of
the following keywords:

 * `AUTHOR` - name of person who created the headerlet
 * `DESCRIP` - short description of the headerlet solution
 * `RMS_RA` - RMS in R.A. at the reference pixel of the WCS stored in the headerlet solution, if updated from the Archive’s default WCS
 * `RMS_DEC` - RMS in Dec. at the reference pixel of the WCS stored in the headerlet solution, if updated from the Archive’s default WCS
 * `NMATCH` - number of sources used in the new solution fit, if updated from the Archive’s default WCS
 * `CATALOG` - astrometric catalog used for headerlet solution
 * `COMMENT` - long description of how the headerlet solution was derived, if updated from Archive’s default WCS

These keywords allow the headerlet to retain enough information about how the
new solution was generated so that a user could determine if it can be applied to his or
her copy of the image.


User-Defined Headerlet
======================
The `headerlet` defined above serves as the default headerlet for any image provided by the HST Archive.  However, should the user perform their own calibrations which they feel improve on the standard calibrations provided by the pipeline, a custom `headerlet` can be provided.  Any `headerlet` should simply include:

    * **Required**: A primary header with specific keywords which specify a unique headerlet name and a targeted image. 
    * **Required**: An SIPWCS extension for each chip which contains the linear WCS as well as any distortion model supported by FITS (for example, updated SIP coefficients)
    * **Optional**: Any additional look up tables with refinements to the polynomial solutions in the SIPWCS extension. Any such extensions should be linked to the SIPWCS extension using the same FITS Distortion Paper conventions used for the NPOLFILE tables. 
    * **Optional**: Detector to image correction array as a separate extension if needed.
    
This custom `headerlet` should be capable of being used to overwrite the existing SCI header WCS keywords to provide a FITS-supported WCS. 


Application of a Headerlet
==========================
Updating an image retrieved from the HST Archive with a `headerlet` only requires a few very simple steps:

    #. Create a headerlet from the original WCS solution in the science image (this step can be turned off).
    #. Delete all WCS information from the science image
    #. Copy the WCS solution from the headerlet to the science observation 
    #. Update the WCSCORR table with the linear WCS keyword values and name of the SIP solution (based on the name of the reference files) from each SIPWCS extension from the `headerlet`, along with the keyword values from the PRIMARY header of the `headerlet`

This process assumes that when an image gets updated with a `headerlet`, the new solution from the `headerlet` should become the prime WCS.  Further implementations of the software to work with `headerlets` can expand on this functionality if necessary.  Initially, the `headerlet` simply needs to be used to update the image's FITS file so that the WCS information can be used at all.

Software Requirements
=====================
Implementing support for the `headerlet` and its use in updating HST FITS files will require a few new software tasks; namely,

- A task which given a science file creates a `headerlet` and writes it to a file.

- A task which given a science file and a `headerlet` applies the `headerlet` to the science file
  
  #. Default behaviour will be to copy the WCS recorded in the `headerlet` as a primary WCS, creating a headerlet with the old solution.

The operation of updating a science file with a `headerlet` only requires the use of basic FITS operations:

- Updating keywords in the science extensions of the file with values from the SIPWCS extensions from the `headerlet`

These operations do not require any computations and can be done using any FITS library. This allows a `headerlet` to be usable by the community even if they do not use the software we develop based on PyFITS and STWCS, both for creating and applying these files.

Headerlet API
=============
This section describes the current draft API for working with `headerlets` as implemented in the `stwcs.wcsutil.headerlet` module.
First, there's a potentially confusing point that should be cleared up:  A `headerlet`, as implemented, is simply a FITS file containing
multiple extensions that contain all the parameters necessary to reproduce the WCS solution in the science image it was created from.
When a `headerlet` is applied to an image, a copy of the original `headerlet` file is appended to the image's HDU list as a special
extension HDU called a `Headerlet HDU`.  A `Headerlet HDU` consists of a simple header describing the `headerlet`, and has as its data
the `headerlet` file itself, (which may be compressed).  A `Headerlet HDU` has an 'XTENSION' value of 'HDRLET'.  Though PyFits can
handle such a non-standard extension type sensibly, this hasn't been tested with other common FITS readers yet.  If it becomes
necessary, `Headerlet HDUs` could be implemented using a standard extension type like 'IMAGE'.

To create a `headerlet` from an image, a `createHeaderlet()` function is provided::

    >>> from stwcs.wcsutil import headerlet
    >>> hdrlet = headerlet.createHeaderlet('j94f05bgq_flt.fits', 'VERSION1')
    >>> type(hdrlet)
    <class 'stwcs.wcsutil.headerlet.Headerlet'>
    >>> hdrlet.info()
    Filename: (No file associated with this HDUList)
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU      12  ()            
    1    SIPWCS      ImageHDU       111  ()            
    2    SIPWCS      ImageHDU       110  ()            
    3    WCSDVARR    ImageHDU        15  (65, 33)      float32
    4    WCSDVARR    ImageHDU        15  (65, 33)      float32
    5    WCSDVARR    ImageHDU        15  (65, 33)      float32
    6    WCSDVARR    ImageHDU        15  (65, 33)      float32
    7    D2IMARR     ImageHDU        12  (4096,)       float32

As you can see, the `Headerlet` object is similar to a normal pyfits `HDUList` object.  `createHeaderlet()` can be given either the path
to a file, or an already open `HDUList` as its first argument.

What do you do with a `Headerlet` object?  Its main purpose is to apply its WCS solution to another file.  This can be done using the
`Headerlet.apply()` method::

    >>> hdrlet.apply('some_other_image.fits')

Or you can use the `applyHeaderlet()` convenience function.  It takes an existing `headerlet` file path or object as its first argument;
the rest of its arguments are the same as `Headerlet.apply()`.  As with `createHeaderlet()` both of these can take a file path or opened
`HDUList` objects as arguments.

When a `headerlet` is applied to an image, an additional `headerlet` containing that image's original WCS solution is automatically created,
and is appended to the file's HDU list as a `Headerlet HDU`.  However, this behavior can be disabled by setting the `createheaderlet` keyword
argument to `False` in either `Headerlet.apply()` or `applyHeaderlet()`.

When opening a file that contains `Headerlet HDU` extensions, it will normally look like this in PyFits::

    >>> import pyfits
    >>> hdul = pyfits.open('94f05bgq_flt_with_hlet.fits')
    >>> hdul.info()
    Filename: j94f05bgq_flt_with_hlet.fits
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU     248  ()            int16
    1    SCI         ImageHDU       286  (4096, 2048)  float32
    2    ERR         ImageHDU        76  (4096, 2048)  float32
    3    DQ          ImageHDU        66  (4096, 2048)  int16
    4    SCI         ImageHDU       282  (4096, 2048)  float32
    5    ERR         ImageHDU        74  (4096, 2048)  float32
    6    DQ          ImageHDU        66  (4096, 2048)  int16
    7    WCSCORR     BinTableHDU     56  10R x 23C     [40A, I, 1A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, D, D, D, D, J, 40A]
    8    WCSDVARR    ImageHDU        15  (65, 33)      float32
    9    WCSDVARR    ImageHDU        15  (65, 33)      float32
    10   WCSDVARR    ImageHDU        15  (65, 33)      float32
    11   WCSDVARR    ImageHDU        15  (65, 33)      float32
    12   D2IMARR     ImageHDU        12  (4096,)       float32
    13   HDRLET  NonstandardExtHDU   13
    14   HDRLET  NonstandardExtHDU   13

The names of the `headerlet` extensions are both HDRLET, but its type shows up as `NonstandardExtHDU`.  Their headers can be read, and while
their data can be read you'd have to know what to do with it (the data is actually either a tar file or a gzipped tar file containing the
`headerlet` file).  However, if you have `stwcs.wcsutil.headerlet` imported, PyFits will recognize these extensions as `Headerlet HDUs`::

    >>> import stwcs.wcsutil.headerlet
    >>> # Note that it's necessary to reopen the file
    >>> hdul = pyfits.open('j94f05bgq_flt_with_hlet.fits')
    >>> hdul.info()
    Filename: j94f05bgq_flt_with_hlet.fits
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU     248  ()            int16
    1    SCI         ImageHDU       286  (4096, 2048)  float32
    2    ERR         ImageHDU        76  (4096, 2048)  float32
    3    DQ          ImageHDU        66  (4096, 2048)  int16
    4    SCI         ImageHDU       282  (4096, 2048)  float32
    5    ERR         ImageHDU        74  (4096, 2048)  float32
    6    DQ          ImageHDU        66  (4096, 2048)  int16
    7    WCSCORR     BinTableHDU     56  10R x 23C     [40A, I, 1A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, D, D, D, D, J, 40A]
    8    WCSDVARR    ImageHDU        15  (65, 33)      float32
    9    WCSDVARR    ImageHDU        15  (65, 33)      float32
    10   WCSDVARR    ImageHDU        15  (65, 33)      float32
    11   WCSDVARR    ImageHDU        15  (65, 33)      float32
    12   D2IMARR     ImageHDU        12  (4096,)       float32
    13   HDRLET      HeaderletHDU    13
    14   HDRLET      HeaderletHDU    13
    >>> print hdul['HDRLET', 1].header.ascard
    XTENSION= 'HDRLET  '           / Headerlet extension                            
    BITPIX  =                    8 / array data type                                
    NAXIS   =                    1 / number of array dimensions                     
    NAXIS1  =               102400 / Axis length                                    
    PCOUNT  =                    0 / number of parameters                           
    GCOUNT  =                    1 / number of groups                               
    EXTNAME = 'HDRLET  '           / name of the headerlet extension                
    HDRNAME = 'j94f05bgq_orig'     / Headerlet name                                 
    DATE    = '2011-04-13T12:14:42' / Date FITS file was generated                  
    SIPNAME = 'IDC_qbu1641sj'      / SIP distortion model name                      
    NPOLFILE= '/grp/hst/acs/lucas/new-npl/qbu16424j_npl.fits' / Non-polynomial corre
    D2IMFILE= '/grp/hst/acs/lucas/new-npl/wfc_ref68col_d2i.fits' / Column correction
    COMPRESS=                    F / Uses gzip compression 

`HeaderletHDU` objects are similar to other HDU objects in PyFits.  However, they have a special `.headerlet` attribute that returns
the actual `headerlet` contained in the HDU data as a `Headerlet` object::

    >>> hdrlet = hdul['HDERLET', 1].headerlet
    >>> hdrlet.info()
    Filename: (No file associated with this HDUList)
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU      12  ()            uint8
    1    SIPWCS      ImageHDU       111  ()            uint8
    2    SIPWCS      ImageHDU       110  ()            uint8
    3    WCSDVARR    ImageHDU        15  (65, 33)      float32
    4    WCSDVARR    ImageHDU        15  (65, 33)      float32
    5    WCSDVARR    ImageHDU        15  (65, 33)      float32
    6    WCSDVARR    ImageHDU        15  (65, 33)      float32
    7    D2IMARR     ImageHDU        12  (4096,)       float32

This is useful if you want to view the contents of the `headerlets` attached to a file.

.. _FITSConventions: http://mediawiki.stsci.edu/mediawiki/index.php/Telescopedia:FITSDistortionConventions
