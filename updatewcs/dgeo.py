from __future__ import division # confidence high

import pyfits
from pytools import fileutil
from stwcs import utils
import numpy as np

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
        
        new_kw = {'DGEOEXT': dgfile}
        return new_kw
    
    updateWCS = classmethod(updateWCS)        

    def applyDgeoCorr(cls, fobj):
        """
        For each science extension in a pyfits file object:
            - create a WCSDVARR extension
            - update science header
            - add/update DGEOEXT keyword
        """
        dgfile = fileutil.osfn(fobj[0].header['DGEOFILE'])
        # Map WCSDVARR EXTVER numbers to extension numbers
        wcsdvarr_ind = cls.getWCSIndex(fobj)
        for ext in fobj:
            try:
                extname = ext.header['EXTNAME'].lower()
            except KeyError:
                continue
            if extname == 'sci':
                extversion = ext.header['EXTVER']
                ccdchip = cls.get_ccdchip(fobj, extname='SCI', extver=extversion)
                binned = utils.getBinning(fobj, extversion)
                header = ext.header
                # get the data arrays from the reference file and transform them for use with SIP
                dx,dy = cls.getData(dgfile, ccdchip)
                idccoeffs = cls.getIDCCoeffs(header)
                if idccoeffs != None:
                    dx, dy = cls.transformData(dx,dy, idccoeffs)
                    
                # Determine EXTVER for the WCSDVARR extension from the DGEO file (EXTNAME, EXTVER) kw.
                # This is used to populate DPj.EXTVER kw
                wcsdvarr_x_version = 2 * extversion -1
                wcsdvarr_y_version = 2 * extversion 
                
                for ename in zip(['DX', 'DY'], [wcsdvarr_x_version,wcsdvarr_y_version],[dx, dy]):
                    cls.addSciExtKw(header, wdvarr_ver=ename[1], dgeo_extname=ename[0])
                    hdu = cls.createDgeoHDU(header, dgeofile=dgfile, \
                        wdvarr_ver=ename[1], dgeo_extname=ename[0], data=ename[2],ccdchip=ccdchip, binned=binned)
                    if wcsdvarr_ind:
                        fobj[wcsdvarr_ind[ename[1]]] = hdu
                    else:
                        fobj.append(hdu)
        
        
    applyDgeoCorr = classmethod(applyDgeoCorr)
              
    def getWCSIndex(cls, fobj):
       
        """
        If fobj has WCSDVARR extensions: 
            returns a mapping of their EXTVER kw to file object extension numbers
        if fobj does not have WCSDVARR extensions:
            an empty dictionary is returned
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
    
    def addSciExtKw(cls, hdr, wdvarr_ver=None, dgeo_extname=None):
        """
        Adds kw to sci extension to define WCSDVARR lookup table extensions
        
        """
        if dgeo_extname =='DX':
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
                
        comments = {cperror: 'Maximum error of dgeo correction for axis %s' % j,  
                    cpdis: 'Prior distortion funcion type',  
                    dpext: 'Version number of WCSDVARR extension containing lookup distortion table', 
                    dpnaxes: 'Number of independent variables in distortion function',
                    dpaxis1: 'Axis number of the jth independent variable in a distortion function', 
                    dpaxis2: 'Axis number of the jth independent variable in a distortion function'
                    }
        
        for key in keys:
            hdr.update(key=key, value=values[key], comment=comments[key], before='HISTORY')
        
    addSciExtKw = classmethod(addSciExtKw)
    
    def getData(cls,dgfile, ccdchip):
        """
        Get the data arrays from the reference DGEO files
        Make sure 'CCDCHIP' in the dgeo file matches "CCDCHIP' in the science file.
        """
        dgf = pyfits.open(dgfile)
        for ext in dgf:
            dgextname  = ext.header.get('EXTNAME',"")
            dgccdchip  = ext.header.get('CCDCHIP',1)
            if dgextname == 'DX' and dgccdchip == ccdchip:
                xdata = ext.data.copy()
                continue
            elif dgextname == 'DY' and dgccdchip == ccdchip:
                ydata = ext.data.copy()
                continue
            else:
                continue
        dgf.close()
        return xdata, ydata
    getData = classmethod(getData)
    
    def transformData(cls, dx, dy, coeffs):
        """
        Transform the DGEO data arrays for use with SIP
        """
        ndx, ndy = np.dot(coeffs, [dx.ravel(), dy.ravel()])
        ndx.shape = dx.shape
        ndy.shape=dy.shape
        return ndx, ndy
    
    transformData = classmethod(transformData)
    
    def getIDCCoeffs(cls, header):
        """
        Return a matrix of the scaled first order IDC coefficients.
        """
        try:
            ocx10 = header['OCX10']
            ocx11 = header['OCX11']
            ocy10 = header['OCY10']
            ocy11 = header['OCY11']
            coeffs = np.array([[ocx11, ocx10], [ocy11,ocy10]], dtype=np.float32)
        except KeyError:
            print 'First order IDCTAB coefficients are not available.\n'
            print 'Cannot convert SIP to IDC coefficients.\n'
            return None
        try:
            idcscale = header['IDCSCALE']
        except KeyError:
            idcscale = 1
            
        return np.linalg.inv(coeffs/idcscale)
    
    getIDCCoeffs = classmethod(getIDCCoeffs)
    
    def createDgeoHDU(cls, sciheader, dgeofile=None, wdvarr_ver=1, dgeo_extname=None,data = None, ccdchip=1, binned=1):
        """
        Creates an HDU to be added to the file object.
        """
        hdr = cls.createDgeoHdr(sciheader, dgeofile=dgeofile, wdvarr_ver=wdvarr_ver, dg_extname=dgeo_extname, ccdchip=ccdchip, binned=binned)
        hdu=pyfits.ImageHDU(header=hdr, data=data)
        return hdu
    
    createDgeoHDU = classmethod(createDgeoHDU)
    
    def createDgeoHdr(cls, sciheader, dgeofile, wdvarr_ver, dg_extname, ccdchip, binned):
        """
        Creates a header for the WCSDVARR extension based on the DGEO reference file 
        and sci extension header. The goal is to always work in image coordinates
        (also for subarrays and binned images. The WCS for the WCSDVARR extension 
        i ssuch that a full size dgeo table is created and then shifted or scaled 
        if the science image is a subarray or binned image.
        """
        dgf = pyfits.open(dgeofile)
        for ext in dgf:
        #for i in range(len(dgf)):
            try:
                dgextname = ext.header['EXTNAME']
                dgextver = ext.header['EXTVER']
            except KeyError:
                continue
            #dgccdchip = ext.header.get('CCDCHIP', 0)
            dgccdchip = cls.get_ccdchip(dgf, extname=dgextname, extver=dgextver)
            if dgextname == dg_extname and dgccdchip == ccdchip:
                dgeo_header = ext.header
                break
            else:
                continue
        dgf.close()
        
        naxis = pyfits.getval(dgeofile, ext=1, key='NAXIS')
        ccdchip = dgextname #dgeo_header['CCDCHIP']
        
        kw = { 'NAXIS': 'Size of the axis', 
                'CRPIX': 'Coordinate system reference pixel', 
                'CRVAL': 'Coordinate system value at reference pixel',
                'CDELT': 'Coordinate increment along axis'}
                
        kw_comm1 = {}
        kw_val1 = {}
        for key in kw.keys():
            for i in range(1, naxis+1):
                si = str(i)
                kw_comm1[key+si] = kw[key]
                
        for i in range(1, naxis+1):
            si = str(i)
            kw_val1['NAXIS'+si] = dgeo_header.get('NAXIS'+si)
            kw_val1['CRPIX'+si] = kw_val1['NAXIS'+si]/2. 
            kw_val1['CDELT'+si] = float(dgeo_header.get('ONAXIS'+si))/ (kw_val1['NAXIS'+si] * binned)
            kw_val1['CRVAL'+si] = (dgeo_header.get('ONAXIS'+si)/2. + \
                                        sciheader.get('LTV'+si, 0.)) / binned
        
                        
        kw_comm0 = {'XTENSION': 'Image extension',
                    'BITPIX': 'IEEE floating point',
                    'NAXIS': 'Number of axes',
                    'EXTNAME': 'WCS distortion array',
                    'EXTVER': 'Distortion array version number',
                    'PCOUNT': 'Special data area of size 0',
                    'GCOUNT': 'One data group',
                    }
        
        kw_val0 = { 'XTENSION': 'IMAGE',
                    'BITPIX': -32,
                    'NAXIS': naxis,
                    'EXTNAME': 'WCSDVARR',
                    'EXTVER':  wdvarr_ver,
                    'PCOUNT': 0,
                    'GCOUNT': 1,
                    'CCDCHIP': ccdchip,
                }
        
        cdl = pyfits.CardList()
        for key in kw_comm0.keys():
            cdl.append(pyfits.Card(key=key, value=kw_val0[key], comment=kw_comm0[key]))
        for key in kw_comm1.keys():
            cdl.append(pyfits.Card(key=key, value=kw_val1[key], comment=kw_comm1[key]))
            
        
        hdr = pyfits.Header(cards=cdl)
        
        return hdr
    
    createDgeoHdr = classmethod(createDgeoHdr)
    
    def get_ccdchip(cls, fobj, extname, extver):
        """
        Given a science file or dgeo file determine CCDCHIP
        """
        ccdchip = 1
        if fobj[0].header['INSTRUME'] == 'ACS' and fobj[0].header['DETECTOR'] == 'WFC':
            ccdchip = fobj[extname, extver].header['CCDCHIP']
        elif fobj[0].header['INSTRUME'] == 'WFC3' and fobj[0].header['DETECTOR'] == 'UVIS':
            ccdchip =  fobj[extname, extver].header['CCDCHIP']
        elif fobj[0].header['INSTRUME'] == 'WFPC2':
            ccdchip = fobj[extname, extver].header['DETECTOR']
        elif fobj[0].header['INSTRUME'] == 'STIS':
            ccdchip = fobj[extname, extver].header['DETECTOR']
        elif fobj[0].header['INSTRUME'] == 'NICMOS':
            ccdchip = fobj[extname, extver].header['CAMERA']
        return ccdchip
        
    get_ccdchip = classmethod(get_ccdchip)
    