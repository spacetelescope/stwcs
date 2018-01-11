.. _appendix2:

Appendix 2 - NPOLFILE Example 
==================================
The NPOLFILE reference file format includes a PRIMARY header describing what kind of 
image should be corrected by this file, along with extensions containing the corrections
for each chip.  

FITS File Extensions
--------------------
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
--------------------

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
----------------------

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


