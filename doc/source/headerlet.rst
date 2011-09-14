.. _headerlet:

**************************************
Headerlet
**************************************
The 'headerlet' serves as a mechanism for encapsulating WCS information
for a single pointing so that it can be used to update the WCS solution of
an image. The concept of a 'headerlet' seeks to provide a solution where
only the WCS solution for an image that has been aligned to an astrometric
catalog can be archived and retrieved for use in updating copies of
that image's WCS information without getting the image data again.
Multiple 'headerlets' could even be provided with each representing the
alignment of an image to a different astrometric solution, giving the end
user the option to get the solution that would allow them to best align
their images with external data of interest to them.  These benefits
can only be realized with the proper definition of a 'headerlet' and
the procedures used to define them and apply them to data.

The headerlet object needs to be as compact as possible while providing
an unambigious and self-consistent WCS solution for an image while
requiring a minimum level of software necessary to apply the headerlet
to an image. 

Headerlet File Structure
-------------------------
This new object complete with the NPOLFILE and the D2IMFILE extensions
derived from the full FITS file fully describes the WCS of each chip
and serves without further modification as the definition of the
`headerlet`. The listing of the FITS extensions for a `headerlet` for
the sample ACS/WFC exposure after writing it out to a file would then be::

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

`DESTIM`   - target image filename (the ROOTNAME keyword of the original archive filename)

`WCSNAME`  - the value of WCSNAME<key> copied from the WCS which was used to create the headerlet

`SIPNAME`  - the name of reference file which contained the original distortion model coefficients. A blank value or 'N/A' will indicate no SIP model was
provided or applied. A value of 'UNKNOWN' indicates a SIP model of unknown origin.

`NPOLFILE` - the name of the NPOLFILE, the reference file which contained the original non-polynomial corrections. The same rules used for SIPNAME apply here as well.

`D2IMFILE` - the name of the D2IMFILE, the reference file which contained the detector to image correction (such as column width correction calibrations). The same rules used for SIPNAME apply here as well.

`DISTNAME` - a concatenation of SIPNAME, NPOLFILE, and D2IMFILE used as a quick reference for the distortion models included with this headerlet.

`UPWCSVER` - version of STWCS used to create the WCS of the original image

`PYWCSVER` - version of PyWCS used to create the WCS of the original image


User-Interface: headerlet
--------------------------
The headerlet module provides those functions necessary 
for creating, updating, and applying headerlets to FITS images.

.. toctree::
   :maxdepth: 2
   
   headerlet_ui
