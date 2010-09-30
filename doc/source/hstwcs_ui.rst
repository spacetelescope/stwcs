.. _hstwcs_ui:

**************************************
HSTWCS UI Examples
**************************************
- Create an HSTWCS object using a pyfits HDUList and an extension number 

  fobj = pyfits.open('some_file.fits')

  w = wcsutil.HSTWCS(fobj, 3)

- Create an HSTWCS object using a qualified file name. 

  w = wcsutil.HSTWCS('j9irw4b1q_flt.fits[sci,1]')

- Create an HSTWCS object using a file name and an extension number. 

  w = wcsutil.HSTWCS('j9irw4b1q_flt.fits', ext=2)
  
- Create an HSTWCS object from WCS with key 'O'.

  w = wcsutil.HSTWCS('j9irw4b1q_flt.fits', ext=2, wcskey='O')

- Create a template HSTWCS object for a DEFAULT object.

  w = wcsutil.HSTWCS(instrument='DEFAULT')