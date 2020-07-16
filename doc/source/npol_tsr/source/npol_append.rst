Appendix
========

Sample WCSDVARR extension header::

    XTENSION= 'IMAGE   '           / Image extension
    BITPIX  =                  -32 / array data type
    NAXIS   =                    2 / number of array dimensions
    NAXIS1  =                   65
    NAXIS2  =                   33
    PCOUNT  =                    0 / number of parameters
    GCOUNT  =                    1 / number of groups
    EXTVER  =                    1 / Distortion array version number
    EXTNAME = 'WCSDVARR'           / WCS distortion array
    CRVAL2  =                  0.0 / Coordinate system value at reference pixel
    CRPIX1  =                  0.0 / Coordinate system reference pixel
    CRPIX2  =                  0.0 / Coordinate system reference pixel
    CRVAL1  =                  0.0 / Coordinate system value at reference pixel
    CDELT1  =                   64 / Coordinate increment along axis
    CDELT2  =                   64 / Coordinate increment along axis


Sample science extension keywords defining a lookup table distortions::

    CPERROR1=                0.058602706   / Maximum error of dgeo correction for axis 1
    CPDIS1  = 'Lookup  '           / Prior distortion function type
    DP1     = 'EXTVER: 1.0'        / Version number of WCSDVARR extension
    DP1     = 'NAXES: 2.0'         / Number of independent variables in CPDIS function
    DP1     = 'AXIS.1: 1.0'        / Axis number of the 1st variable in a CPDIS function
    DP1     = 'AXIS.2: 2.0'        / Axis number of the 2nd variable in a CPDIS function
    CPERROR2=        0.072911568  / Maximum error of dgeo correction for axis 2
    CPDIS2  = 'Lookup  '           / Prior distortion function type
    DP2     = 'EXTVER: 2.0'        / Version number of WCSDVARR extension
    DP2     = 'NAXES: 2.0'         / Number of independent variables in CPDIS function
    DP2     = 'AXIS.1: 1.0'        / Axis number of the 1st variable in a CPDIS function
    DP2     = 'AXIS.2: 2.0'        / Axis number of the 2nd variable in a CPDIS function
