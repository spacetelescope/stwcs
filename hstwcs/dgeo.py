import pyfits
from pytools import fileutil
#from hstwcs.mappings import dgeo_vals
import numpy

class DGEOCorr(object):
    """
    Purpose
    =======
    Defines a Lookup table prior distortion correction as per WCS paper IV.
    It uses a reference file defined by the DGEOFILE keyword in the primary header.
    
    Algorithm
    =========
    - Using extensions in the reference file create a WCSDVARR extension 
      and add it to the file object.
    - Add record-valued keywords which describe the lookup tables to the 
      science extension header
    - Add a keyword 'DGEOFILE' to the science extension header, whose
      value is the reference file used to create the WCSVARR extension
    
    If WCSDVARR extensions exist, subsequent updates will overwrite them. 
    If not, they will be added to the file object.
    
    It is assumed that the DGEO reference files were created to work with IDC tables
    but will be applied with SIP coefficients. A transformation is applied to correct 
    for the fact that the lookup tables will be applied before the first order coefficients
    which are in the CD matrix when the SIP convention is used.
    """
    
    def updateWCS(cls, fobj):
        """
        :Parameters:
        `fobj`: pyfits object
                Science file, for which a distortion correction in a DGEOFILE is available
                
        """
        assert isinstance(fobj, pyfits.NP_pyfits.HDUList)
        cls.applyDgeoCorr(fobj)
        dgfile = fobj[0].header['DGEOFILE']
        
        new_kw = {'DGEOFILE': dgfile}
        return new_kw
    
    updateWCS = classmethod(updateWCS)        

    def applyDgeoCorr(cls, fobj):
        """
        For each science extension in a pyfits file object:
            - create a WCSDVARR extension
            - update science header
            - add/update DGEOFILE keyword
        """
        dgfile = fileutil.osfn(fobj[0].header['DGEOFILE'])
        instrument = fobj[0].header.get('INSTRUME', None)
        # Map WCSDVARR EXTVER numbers to extension numbers
        wcsdvarr_ind = cls.getWCSIndex(fobj)
        for ext in fobj:
            try:
                extname = ext.header['EXTNAME'].lower()
            except KeyError:
                continue
            if extname == 'sci':
                extversion = ext.header['EXTVER']
                header = ext.header
                # get the data arrays from the reference file and transform them for use with SIP
                dx,dy = cls.getData(dgfile, extversion)
                ndx, ndy = cls.transformData(header, dx,dy)
                # determine EXTVER for the WCSDVARR extension from the DGEO file (EXTNAME, EXTVER) kw
                wcsdvarr_x_version = 2 * extversion -1
                wcsdvarr_y_version = 2 * extversion 
                
                for ename in zip(['DX', 'DY'], [wcsdvarr_x_version,wcsdvarr_y_version],[ndx, ndy]):
                    cls.addSciExtKw(header, wdvarr_ver=ename[1], dgeo_name=ename[0])
                    hdu = cls.createDgeoHDU(header, dgeofile=dgfile, wdvarr_ver=ename[1],dgeo_name=ename[0], data=ename[2],extver=extversion)
                    if wcsdvarr_ind:
                        fobj[wcsdvarr_ind[ename[1]]] = hdu
                    else:
                        fobj.append(hdu)
        
        
    applyDgeoCorr = classmethod(applyDgeoCorr)
              
    def getWCSIndex(cls, fobj):
        """
        If fobj has WCSDVARR extensions: 
            returns a mapping of their EXTVER kw are mapped to extension numbers
        if fobj does not have WCSDVARR extensions:
            an empty dictionary is returned.
        """
        wcsd = {}
        for e in range(len(fobj)):
            try:
                ename = fobj[e].header['EXTNAME']
            except KeyError:
                continue
            if ename == 'WCSDVARR':
                wcsd[fobj[e].header['EXTVER']] = e
        
        return wcsd
        
    getWCSIndex = classmethod(getWCSIndex)
    
    def addSciExtKw(cls, hdr, wdvarr_ver=None, dgeo_name=None):
        """
        Adds kw to sci extension to define WCSDVARR lookup table extensions
        
        """
        if dgeo_name =='DX':
            j=1
        else:
            j=2
        
        cperror = 'CPERROR%s' %j
        cpdis = 'CPDIS%s' %j
        dpext = 'DP%s.' %j + 'EXTVER'
        dpnaxes = 'DP%s.' %j +'NAXES'
        dpaxis1 = 'DP%s.' %j+'AXIS.1'
        dpaxis2 = 'DP%s.' %j+'AXIS.2'
        keys = [cperror, cpdis, dpext, dpnaxes, dpaxis1, dpaxis2]
        values = {cperror: 0.0, cpdis: 'Lookup',  dpext: wdvarr_ver, dpnaxes: 2,
                dpaxis1: 1, dpaxis2: 2}
                
        comments = {cperror: 'Maximum error of dgeo correction for axis %s' % (wdvarr_ver/2), 
                    cpdis: 'Prior distortion funcion type',  
                    dpext: 'Version number of WCSDVARR extension containing lookup distortion table', 
                    dpnaxes: 'Number of independent variables in distortion function',
                    dpaxis1: 'Axis number of the jth independent variable in a distortion function', 
                    dpaxis2: 'Axis number of the jth independent variable in a distortion function'
                    }
        
        for key in keys:
            hdr.update(key=key, value=values[key], comment=comments[key], before='HISTORY')
        
    addSciExtKw = classmethod(addSciExtKw)
    
    def getData(cls,dgfile, extver):
        """
        Get the data arrays from the reference DGEO files
        """
        xdata = pyfits.getdata(dgfile, ext=('DX',extver))
        ydata = pyfits.getdata(dgfile, ext=('DY',extver))
        return xdata, ydata
    getData = classmethod(getData)
    
    def transformData(cls, header, dx, dy):
        """
        Transform the DGEO data arrays for use with SIP
        """
        coeffs = cls.getCoeffs(header)
        idcscale = header['IDCSCALE']
        sclcoeffs = numpy.linalg.inv(coeffs/idcscale)
        ndx, ndy = numpy.dot(sclcoeffs, [dx.ravel(), dy.ravel()])
        ndx.shape = dx.shape
        ndy.shape=dy.shape
        return ndx, ndy
    transformData = classmethod(transformData)
    
    def getCoeffs(cls, header):
        """
        Return a matrix of the first order IDC coefficients.
        """
        try:
            ocx10 = header['OCX10']
            ocx11 = header['OCX11']
            ocy10 = header['OCY10']
            ocy11 = header['OCY11']
        except AttributeError:
            print 'First order IDCTAB coefficients are not available.\n'
            print 'Cannot convert SIP to IDC coefficients.\n'
            return None
        return numpy.array([[ocx11, ocx10], [ocy11,ocy10]], dtype=numpy.float32)
    
    getCoeffs = classmethod(getCoeffs)
    
    def createDgeoHDU(cls, sciheader, dgeofile=None, wdvarr_ver=1, dgeo_name=None,data = None, extver=1):
        """
        Creates an HDU to be added to the file object.
        """
        hdr = cls.createDgeoHdr(sciheader, dgeofile=dgeofile, wdvarr_ver=wdvarr_ver, dgeoname=dgeo_name, extver=extver)
        hdu=pyfits.ImageHDU(header=hdr, data=data)
        return hdu
    
    createDgeoHDU = classmethod(createDgeoHDU)
    
    def createDgeoHdr(cls, sciheader, dgeofile, wdvarr_ver, dgeoname, extver):
        """
        Creates a header for the WCSDVARR extension based on the DGEO reference file 
        and sci extension header.
        """
        dgeo_header = pyfits.getheader(dgeofile, ext=(dgeoname,extver))
        sci_naxis1 = sciheader['NAXIS1']
        sci_naxis2 = sciheader['NAXIS2']
        sci_crpix1 = sciheader['CRPIX1']
        sci_crpix2 = sciheader['CRPIX2']
        
        naxis1 = dgeo_header['naxis1']
        naxis2 = dgeo_header['naxis2'] 
        extver = dgeo_header['extver'] 
        crpix1 = naxis1/2.
        crpix2 = naxis2/2.
        cdelt1 = sci_naxis1/naxis1 
        cdelt2 = sci_naxis2/naxis2 
        crval1 = sci_crpix1
        crval2 = sci_crpix2
        keys = ['XTENSION','BITPIX','NAXIS','NAXIS1','NAXIS2',
              'EXTNAME','EXTVER','PCOUNT','GCOUNT','CRPIX1',
                        'CDELT1','CRVAL1','CRPIX2','CDELT2','CRVAL2']
                        
        comments = {'XTENSION': 'Image extension',
                    'BITPIX': 'IEEE floating point',
                    'NAXIS': 'Number of image columns',
                    'NAXIS1': 'Number of image columns',
                    'NAXIS2': 'Number of image rows',
                    'EXTNAME': 'WCS distortion array',
                    'EXTVER': 'Distortion array version number',
                    'PCOUNT': 'Special data area of size 0',
                    'GCOUNT': 'One data group',
                    'CRPIX1': 'Distortion array reference pixel',
                    'CDELT1': 'Grid step size in first coordinate',
                    'CRVAL1': 'Image array pixel coordinate',
                    'CRPIX2': 'Distortion array reference pixel',
                    'CDELT2': 'Grid step size in second coordinate',
                    'CRVAL2': 'Image array pixel coordinate'}
        
        values = {'XTENSION': 'IMAGE',
                'BITPIX': -32,
                'NAXIS': 2,
                'NAXIS1': naxis1,
                'NAXIS2': naxis2,
                'EXTNAME': 'WCSDVARR',
                'EXTVER':  wdvarr_ver,
                'PCOUNT': 0,
                'GCOUNT': 1,
                'CRPIX1': crpix1,
                'CDELT1': cdelt1,
                'CRVAL1': crval1,
                'CRPIX2': crpix1,
                'CDELT2': cdelt2,
                'CRVAL2': crval2
                }
                    
        
        cdl = pyfits.CardList()
        for c in keys:
            cdl.append(pyfits.Card(key=c, value=values[c], comment=comments[c]))

        hdr = pyfits.Header(cards=cdl)
        return hdr
    
    createDgeoHdr = classmethod(createDgeoHdr)
    
