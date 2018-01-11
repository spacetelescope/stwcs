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
A full example of the keywords added to a ``SCI`` extension header is presented in :ref:`appendix1`.

A complete description of the conversion of the DGEOFILE reference data into NPOLFILE reference
files can be found in the report on the ``npolfile-tsr``.


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
reference file also remains very simple, as illustrated in :ref:`appendix2`.

Applying these corrections starts by reading the two 65 x 33 
arrays into memory with each input ACS/WFC chip WCS (one for 
X offsets and one for Y offsets). Bi-linear interpolation based on the input pixel 
position then gets used on-the-fly to extract the final offset from this reference 
file. Initial versions of these sub-sampled NPOLFILE reference files for ACS have 
been derived from the current full-size DGEOFILEs, and testing indicates residuals 
only on the order of 0.02 pixels or less remain when compared to the original calibration. 

