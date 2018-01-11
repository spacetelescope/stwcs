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
:ref:`appendix3`. 

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

