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


