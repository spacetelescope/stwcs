===================================================================
Definition of a Headerlet and It's Role in Updating WCS Information
===================================================================
:Authors: Warren Hack, Nadezhda Dencheva

Abstract
========
The 'headerlet' serves as a mechanism for encapsulating WCS information which can be used to update the WCS solution of an image. This object needs to be as compact as possible while providing an unambigious and self-consistent WCS solution for an image while requiring a minimum level of software necessary to apply the headerlet to an image.  These basic requirements come from the desire to create a mechanism for passing improved astrometric solutions for HST data and provide those solutions in a manner that would not require getting entirely new images from an archive when only the WCS information has been updated. This report describes the format and contents of a headerlet along with the procedures which would be used to update images with the updated WCS information from the headerlet.   

Introduction
============
The HST Archive provides access to all HST data, most of which can not be readily combined together due to errors in guide star astrometry imposing offsets between images taken using different pairs of guide stars.  A lot of effort has gone into computing those offsets so that all the images taken at a particular pointing can be combined successfully with astrometry matched to external astrometric catalogs. Unfortunately, there is no current mechanism for passing those updated solutions along to the community without providing entirely new copies of all the data.  

The concept of a 'headerlet' seeks to provide a solution where only the WCS solution for an image that has been aligned to an astrometric catalog can be archived and retrieved for use in updating copies of that image's WCS information without getting the image data again.  Multiple 'headerlets' could even be provided with each representing the alignment of an image to a different astrometric solution, giving the end user the option to get the solution that would allow them to best align their images with external data of interest to them.  These benefits can only be realized with the proper definition of a 'headerlet' and the procedures used to define them and apply them to data. 

Source Image
============
All users get their copy of an image from the HST Archive (OTFR) after being processed with the latest image calibrations, including applying the latest available distortion models. The calibrated image gets passed to the user with the ``_flt.fits`` suffix and is referred to as the ``FLT`` image.  These FLT images serve as the inputs to MultiDrizzle in order to apply the distortion models and combine the images into a single ``DRZ`` product.  

The WCS information in the FLT images has been updated to include the full distortion model, including the full polynomial solution from the IDCTAB and all the corrections formerly combined into the DGEOFILE. The FITSConventions_ report by Dencheva (currently available online) contains the full description of the conventions used to describe all these components in a FITS file. The header now contains the following set of keywords and extensions to fully describe the WCS with distortion:

  * **Linear WCS keywords**: specifically, CRPIX, CRVAL, CTYPE, CD matrix keywords
  * **SIP coefficients**: A_*_* and B_*_*, A_ORDER, B_ORDER, 
    OCX10, OCX11, OCY10, and OCY11 keywords
  * **NPOL file**: if an NPOLFILE has been specified for the image, 
    CPDIS and DP keywords to point to WCSDVARR extensions (Paper IV convention)
  * **Column correction file**: if a D2IMFILE has been specified for use with the image, 
    the D2IMEXT and D2IMERR keywords signify the use of a D2IM file extension
  * **WCSDVARR extensions**: 2 extensions for each chip with lookup tables containing 
    the NPOLFILE corrections, with each extension corresponding to an axis of 
    the image (X correction or Y correction) (if an NPOLFILE has been applied to the image)
  * **D2IMVARR extension**: an extension with a lookup table containing the 
    column-correction for that chip from the D2IMFILE.
 

Each science header will have its own set of these keywords and extensions that will be kept together as part of the headerlet definition.  This avoids any ambiguity as to what solution was used for any given WCS. A summary of all the WCS solutions for all the chips can be written out to the WCSCORR table.

An ACS/WFC exposure would end up with the following set of extensions::

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
    12      BINTABLE    WCSCORR                     18Fx10R

There may be a lot of extensions appended to this FITS file, but the sum total of all these new extensions comes to approximately 100kB for ACS/WFC images, making them a fairly space efficient means of managing all the distortion and WCS information. 

Headerlet Definition
====================
The `headerlet` needs to be a self-consistent, fully described definition of a WCS and its distortion.  The WCS and SIP coefficients get derived from the SCI header directly, along with the keywords which refer to the extensions with the NPOLFILE and D2IMFILE corrections.  This original WCS information in general will not be deleted from the file so that the user can always revert back to the original WCS solution at any time. 
There will be a way to delete it in cases when the user want to update the model and does not need to keep the old model.

New WCS Extension
-----------------
A new extension, named SIPWCS, containing all the WCS-related keywords from the SCI header, including all keywords referring to NPOL and D2IM extensions as well as all sets of alternate WCS keywords, can be created to serve as the record of the original WCS. Keywords (TBD) recording the alignment information are recorded in this header as well. All the sets of linear WCS keywords stored using FITS Paper I Multiple WCS Standard would be defined using the same set of distortion coefficients written to the SIP keywords and NPOL files.  This insures that all the information in the header remains consistent. The keywords in this extension can be used to overwrite the keywords in the corresponding SCI header to update the WCS solution for each chip without any further modification or computation. The new extension then serves not only as a record of all the WCS solutions derived for the image, but also the source of values for restoring the SCI header WCS when desired.  


Updated Example FITS Listing
-----------------------------
The updated listing of the sample ACS/WFC exposure FITS file extensions would then be::

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
    12      IMAGE       SIPWCS                1                  8
    13      IMAGE       SIPWCS                2                  8
    14      BINTABLE    WCSCORR                     18Fx10R

This new extension along with the NPOLFILE and the D2IMFILE extensions fully describe the WCS of each chip and can serve without further modification as the definition of the `headerlet`. The listing of the FITS extensions for a `headerlet` for the sample ACS/WFC exposure after writing it out to a file would then be::

    EXT#  FITSNAME      FILENAME              EXTVE DIMENS       BITPI OBJECT       

    0     j8hw27c4q     j8hw27c4q_hdr.fits                       16
    1       IMAGE       D2IMARR               1     4096         -32                
    2       IMAGE       WCSDVARR              1     64x32        -32                
    3       IMAGE       WCSDVARR              2     64x32        -32                
    4       IMAGE       WCSDVARR              3     64x32        -32                
    5       IMAGE       WCSDVARR              4     64x32        -32                
    6       IMAGE       SIPWCS                1                  8
    7       IMAGE       SIPWCS                2                  8

This file now fully describes the WCS solution for this image, complete with all the distortion information used to originally define the solution. No further reference files or computations would be needed when this `headerlet` gets used to update an image.

The primary header must have 4 required keywords:

`HDRNAME`  - a unique name for the headerlet
`DISTIM`   - target image filename (the original archive filename)
`STWCSVER` - version of STWCS used to create the headerlet
`PYWCSVER` - version of PyWCS used to create the headerlet

User-Defined Headerlet
======================
The `headerlet` defined above serves as the default headerlet for any image provided by the HST Archive.  However, should the user perform their own calibrations which they feel improve on the standard calibrations provided by the pipeline, a custom `headerlet` can be provided.  Any `headerlet` should simply include:

    * **Required**: A primary header with specific keywords which specify a unique headerlet name and a targeted image. 
    * **Required**: An SIPWCS extension for each chip which contains the linear WCS as well as any distortion model supported by FITS (for example, updated SIP coefficients)
    * **Optional**: Any additional look up tables with refinements to the polynomial solutions in the SIPWCS extension. Any such extensions should be linked to the SIPWCS extension using the same Paper IV conventions used for the NPOLFILE tables. 
    * **Optional**: Detector to image correction array as a separate extension if needed.
    
This custom `headerlet` should be capable of being used to overwrite the existing SCI header WCS keywords to provide a FITS-supported WCS. 


Application of a Headerlet
==========================
Updating an image retrieved from the HST Archive with a `headerlet` only requires a few very simple steps:

    #. Append the `headerlet` to the FITS file
    #. Update the extver IDs for the NPOLFILE and D2IMFILE keywords in the headerlet SIPWCS extensions to point to the actual extver values for the new extensions
    #. Overwrite the SCI header keywords for each chip with the same keywords from the SIPWCS extension that corresponds to the same chip from the newly appended `headerlet`
    #. Add a keyword `SIPVER` to each science header with a value of the appropriate SIPWCS' `EXTVER` keyword.
    #. Update the WCSCORR table with the linear WCS keyword values and name of the SIP solution from each SIPWCS extension from the `headerlet`

This process assumes that when an image gets updated with a `headerlet`, the new solution from the `headerlet` should become the prime WCS.  Further implementations of the software to work with `headerlets` can expand on this functionality if necessary.  Initially, the `headerlet` simply needs to be used to update the image's FITS file so that the WCS information can be used at all.

Software Requirements
=====================

- A task which given a science file creates a headerlet and writes it to a file.

- A task which given a science file and a headerlet applies the headerlet to the science file
  
  #. Default behaviour will be to append the headerlet to the file and copy the WCS recorded in the headerlet as a primary WCS.
  #. It will be possible (optionally) to copy the updated science file to a new file and keep the original science file locally unchanged.

.. _FITSConventions: http://mediawiki.stsci.edu/mediawiki/index.php/Telescopedia:FITSDistortionConventions
