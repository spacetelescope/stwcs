.. _appendix3:


Appendix 3 - D2IMFILE Example 
==================================

The D2IMFILE reference file only contains a single 1-D array that should correct the
column (row) values based on the value of the 'AXISCORR' keyword in the SCI header. 

FITS File Extensions
--------------------

This simple reference file, therefore, contains only 2 extensions; namely,

::

 Filename: /grp/hst/cdbs/jref/v971826mj_d2i.fits
 No.    Name         Type      Cards   Dimensions   Format
 0    PRIMARY     PrimaryHDU      35   ()           int16   
 1    DX          ImageHDU        18   (4096,)      float32   

PRIMARY Header
--------------------

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
----------------------

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


.. [DistortionPaper] Calabretta M. R., Valdes F. G., Greisen E. W., and Allen S. L., 2004, 
    "Representations of distortions in FITS world coordinate systems",[cited 2012 Sept 18], 
    Available from: http://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf

.. [SIPConvention] Shupe D.L., Hook R.N., 2008, "The SIP Convention for Representing Distortion in FITS Image
    Headers", [cited 2012 Sept 18], Available from: http://fits.gsfc.nasa.gov/registry/sip.html

