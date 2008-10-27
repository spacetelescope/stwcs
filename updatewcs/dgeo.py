import pyfits
from pytools import fileutil
from hstwcs.mappings import dgeo_vals

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
    - Add record-valued keywords which describe the looikup tables to the 
      science extension header
    - Add a keyword 'DGEOFILE' to the science extension header, whose
      value is the reference file used to create the WCSVARR extension
    
    If WCSDVARR extensions exist, subsequent updates will overwrite them. 
    If not, they will be added to the file object.
    
    """
    def updateWCS(cls, fobj):
        """
        :Parameters:
        `fobj`: pyfits object
                Science file, for which a distortion correc tion in a DGEOFILE is available
                
        """
        assert isinstance(fobj, pyfits.NP_pyfits.HDUList)
        #self.fobj = fobj
        cls.applyDgeoCorr(fobj)
        dgfile = fobj[0].header['DGEOFILE']
        
        new_kw = {'DGEOFIEL': dgfile}
        return new_kw
    
    updateWCS = classmethod(updateWCS)        

    def applyDgeoCorr(cls,fobj):
        """
        For each science extension in a pyfits file object:
            - create a WCSDVARR extension
            - update science header
            - add/update DGEOFILE keyword
        """
        
        dgeover = 0
        dgfile = fileutil.osfn(fobj[0].header['DGEOFILE'])
        instrument = fobj[0].header.get('INSTRUME', None)
        wcsdvarr_ind = cls.getWCSIndex(fobj)
        for ext in fobj:
            try:
                extname = ext.header['EXTNAME'].lower()
            except KeyError:
                continue
            if extname == 'sci':
                extversion = ext.header['EXTVER']
                for ename in ['DX', 'DY']:
                    dgeover +=1
                              
                    cls.addSciExtKw(ext.header, extver=dgeover, ename=ename)
                    hdu = cls.createDgeoHDU(instrument, dgeofile=dgfile, dgeover=dgeover,ename=ename, extver=extversion)
                    if wcsdvarr_ind:
                        fobj[wcsdvarr_ind[dgeover]] = hdu
                    else:
                        fobj.append(hdu)
                
    applyDgeoCorr = classmethod(applyDgeoCorr)
                
    def getWCSIndex(cls, fobj):
        """
        Returns a mapping of WCSDVARR 'EXTVER' and pyfits.HDUList indeces.
        
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
    
    def addSciExtKw(cls, hdr, extver=None, ename=None):
        """
        Adds kw to sci extension to define dgeo correction extensions
        kw to be added to sci ext:
        CPERRORj
        CPDISj
        DPj.EXTVER
        DPj.NAXES
        DPj.AXIS.i (i=DP1.NAXES)
        """
        if ename =='DX':
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
        values = {cperror: 0.0, cpdis: 'Lookup',  dpext: extver, dpnaxes: 2,
                dpaxis1: 1, dpaxis2: 2}
                
        comments = {cperror: 'Maximum error of dgeo correction for axis %s' % (extver/2), 
                    cpdis: 'Prior distortion funcion type',  
                    dpext: 'Version number of WCSDVARR extension containing lookup distortion table', 
                    dpnaxes: 'Number of independent variables in distortion function',
                    dpaxis1: 'Axis number of the jth independent variable in a distortion function', 
                    dpaxis2: 'Axis number of the jth independent variable in a distortion function'
                    }
        
        for key in keys:
            hdr.update(key=key, value=values[key], comment=comments[key], before='HISTORY')
        
    addSciExtKw = classmethod(addSciExtKw)
    
    def getDgeoData(cls, dgfile=None, ename=None, extver=1):
        """
        Given a dgeo file name, creates an array to be used 
        as a data array in the dgeo extension.
        """ 
        return pyfits.getdata(dgfile, ext=(ename,extver))

    getDgeoData = classmethod(getDgeoData)
            
    def createDgeoHDU(cls, instrument, dgeofile=None, dgeover=1, ename=None,extver=1):
        """
        Creates an HDU to be added to the file object.
        """
        dgeokw = {'naxis1':None, 'naxis2':None, 'extver':dgeover, 'crpix1':None, 
                'crpix2':None, 'cdelt1':None, 'cdelt2':None, 'crval1':None, 'crval2':None}
        hdr = cls.createDgeoHdr(instrument, **dgeokw)
        data = cls.getDgeoData(dgfile=dgeofile, ename=ename, extver=extver)
        hdu=pyfits.ImageHDU(header=hdr, data=data)
        return hdu
    createDgeoHDU = classmethod(createDgeoHDU)
    
    def createDgeoHdr(cls, instr, **kw):
        """
        Creates a header for the dgeo extension based on dgeo file 
        and sci extension header.
        **kw = {'naxis1':None, 'naxis2':None, 'extver':None, 'crpix1':None, 
                    'crpix2':None, 'cdelt1':None, 'cdelt2':None, 'crval1':None, 'crval2':None}
        """
        #instr = self.fobj[0].header['INSTRUME']
        instr_vals = dgeo_vals[instr]
        naxis1 = kw['naxis1'] or instr_vals['naxis1']
        naxis2 = kw['naxis2'] or instr_vals['naxis2']
        extver = kw['extver'] or instr_vals['extver']
        crpix1 = kw['crpix1'] or instr_vals['crpix1']
        crpix2 = kw['crpix2'] or instr_vals['crpix2']
        cdelt1 = kw['cdelt1'] or instr_vals['cdelt1']
        cdelt2 = kw['cdelt2'] or instr_vals['cdelt2']
        crval1 = kw['crval1'] or instr_vals['crval1']
        crval2 = kw['crval2'] or instr_vals['crval2']
            
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
                'EXTVER':  extver,
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
    
