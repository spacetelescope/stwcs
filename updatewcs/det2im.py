from __future__ import division # confidence high

import pyfits
from pytools import fileutil
import utils

class DET2IMCorr(object):
    def updateWCS(cls, fobj):
        """
        Parameters
        ----------
        fobj: pyfits object
              Science file, for which a detector to image correction 
              is available
              
        Notes
        -----
        Uses the file pointed to in the primary header keyword 'D2IMFILE' 
        to create an extension with a detector to image correction. 
        """
        assert isinstance(fobj, pyfits.HDUList)
        
        d2imfile = fileutil.osfn(fobj[0].header['D2IMFILE'])
        axiscorr = cls.getAxisCorr(d2imfile)
        d2imerr = pyfits.getdata(d2imfile, ext=1).max()
        if axiscorr == None:
            new_kw = {}
        else:
            new_kw = {'D2IMEXT': d2imfile, 'AXISCORR': axiscorr, 'D2IMERR': d2imerr}
            cls.applyDet2ImCorr(fobj,axiscorr)
        cls.updatehdr(fobj, new_kw)
    
    updateWCS = classmethod(updateWCS)        
    
    def getAxisCorr(cls, refname):
        try:
            direction = pyfits.getval(refname, ext=1, key='EXTNAME')
            if direction == 'DX': return 1
            elif direction == 'DY': return 2
            else: 
                print '\tDET2IM correction expects the reference file to have'
                print '\tan EXTNAME keyword of value "DX"  or "DY"'
                return None
        except AttributeError:
            print "\tAxis to which to apply the detector to image "
            print "\tcorrection cannot be determined because the reference "
            print "\tfile %s is missing a keyword EXTNAME" % refname
            direction = None
        return direction
    getAxisCorr = classmethod(getAxisCorr)
    
    def applyDet2ImCorr(cls,fobj, axiscorr):
        binned = utils.getBinning(fobj)
        hdu = cls.createDgeoHDU(fobj, axiscorr, binned)
        d2imarr_ind = cls.getD2imIndex(fobj)
        if d2imarr_ind:
            fobj[d2imarr_ind] = hdu
        else:
            fobj.append(hdu)
    applyDet2ImCorr = classmethod(applyDet2ImCorr)
    
    def getD2imIndex(cls,fobj):
        index = None
        for e in range(len(fobj)):
            try:
                ename = fobj[e].header['EXTNAME']
            except KeyError:
                continue
            if ename == 'D2IMARR':
                index = e
        return index
    getD2imIndex = classmethod(getD2imIndex)
    
    def createDgeoHDU(cls, fobj, axiscorr, binned=1):
        d2imfile = fileutil.osfn(fobj[0].header['D2IMFILE'])
        d2im_data = pyfits.getdata(d2imfile, ext=1)
        sci_hdr = fobj['sci',1].header
        d2im_hdr = cls.createDet2ImHdr(fobj, binned)
        hdu = pyfits.ImageHDU(header=d2im_hdr, data=d2im_data)
        
        return hdu
    
    createDgeoHDU = classmethod(createDgeoHDU)
    
    def createDet2ImHdr(cls, fobj, binned=1):        
        """
        Creates a header for the D2IMARR extension based on the 
        reference file recorded in D2IMFILE keyword in the primary header.
        fobj - the science  file
        
        """
        d2imfile = fileutil.osfn(fobj[0].header['D2IMFILE'])
        axiscorr = cls.getAxisCorr(d2imfile)
        sci_hdr = fobj[1].header
        data_shape = pyfits.getdata(d2imfile, ext=1).shape
        naxis = pyfits.getval(d2imfile, ext=1, key='NAXIS')
        
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
            kw_val1['NAXIS'+si] = data_shape[i-1]
            kw_val1['CRPIX'+si] = data_shape[i-1]/2.
            kw_val1['CDELT'+si] = 1./binned
            kw_val1['CRVAL'+si] = (sci_hdr.get('NAXIS'+si, 1)/2. + \
                                        sci_hdr.get('LTV'+si, 0.)) / binned
        
                        
        kw_comm0 = {'XTENSION': 'Image extension',
                    'BITPIX': 'IEEE floating point',
                    'NAXIS': 'Number of axes',
                    'EXTNAME': 'WCS distortion array',
                    'EXTVER': 'Distortion array version number',
                    'PCOUNT': 'Special data area of size 0',
                    'GCOUNT': 'One data group',
                    'AXISCORR': 'Direction in which the det2im correction is applied'}
        
        kw_val0 = { 'XTENSION': 'IMAGE',
                    'BITPIX': -32,
                    'NAXIS': naxis,
                    'EXTNAME': 'D2IMARR',
                    'EXTVER':  1,
                    'PCOUNT': 0,
                    'GCOUNT': 1,
                    'AXISCORR': axiscorr
                }
                    
        
        cdl = pyfits.CardList()
        for key in kw_comm0.keys():
            cdl.append(pyfits.Card(key=key, value=kw_val0[key], comment=kw_comm0[key]))
        for key in kw_comm1.keys():
            cdl.append(pyfits.Card(key=key, value=kw_val1[key], comment=kw_comm1[key]))
            
        hdr = pyfits.Header(cards=cdl)
        return hdr
    
    createDet2ImHdr = classmethod(createDet2ImHdr)
    
    def updatehdr(cls, fobj, kwdict):
        """
        Update extension headers to keep record of the files used for the 
        detector to image correction.
        """
        for ext in fobj:
            try:
                extname = ext.header['EXTNAME'].lower()
            except KeyError:
                continue
            if extname == 'sci':
                for kw in kwdict:
                    ext.header.update(kw, kwdict[kw])
            else:
                continue
    updatehdr = classmethod(updatehdr)
    
