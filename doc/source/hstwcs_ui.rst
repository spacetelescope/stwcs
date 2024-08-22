.. _hstwcs_ui:

**************************************
HSTWCS Examples
**************************************

--------------------------
Create an HSTWCS Object
--------------------------
- Create an HSTWCS object using a pyfits HDUList and an extension number 

  ``from astropy.io import fits``

  ``from stwcs import wcsutil``

  ``hdu = fits.open('some_file.fits')``

  ``w = wcsutil.HSTWCS(hdu, 3)``

- Create an HSTWCS object using a qualified file name. 

  ``w = wcsutil.HSTWCS('j9irw4b1q_flt.fits[sci,1]')``

- Create an HSTWCS object using a file name and an extension number. 

  ``w = wcsutil.HSTWCS('j9irw4b1q_flt.fits', ext=2)``
  
- Create an HSTWCS object from WCS with key 'O'.

  ``w = wcsutil.HSTWCS('j9irw4b1q_flt.fits', ext=2, wcskey='O')``
  
----------------------------------
Coordinate Transformation Examples
----------------------------------
All coordinate transformation functions accept input coordinates 
as 2D numpy arrays or 2 sequences of X and Y coordinates. 

``inpix = np.array([[1., 2.], [1,3], [1,4], [1,5]])``

or

``X = [1.,1.,1.,1.]``

``Y = np.array([2.,3.,4.,5.])``

In addition all transformation functions require an `origin` parameter 
which specifies if the coordinates are 0 or 1 based. For example in FITS 
and Fortran, coordinates start from 1, while in Python and C, the index 
of the first image pixel is (0,0).

- Apply the entire detector to sky transformation at once:

 ``outpix = w.all_pix2world(inpix, 1)``

 ``outpix = w.all_pix2world(X, Y, 1)``

- The same transformation can be done in separate steps:

1. Apply the detector to image correction

 ``dpx = w.det2im(inpix,1)``

2. Aply the SIP polynomial distortion

 ``spx = w.sip_pix2foc(dpx, 1)``
 
3. Apply the non-polynomial distortion from the lookup table

 ``lutpx = w.p4_pix2foc(dpx,1)``
 
4. The undistorted coordinates are the sum of the input coordinates with 
   the deltas for the distortion corrections.

 ``fpix = dpx + (spx-dpx) +(lutpx-dpx)``
 
5. Finally the transformation from undistorted to world coordinates is done 
   by applying the linear WCS.
 
 ``wpix = w.wcs_pix2world(fpix, 1)``
 
 
