import pyfits
from pytools import fileutil

class DET2IMCorr(object):
    def updateWCS(cls, fobj):
        """
        :Parameters:
        `fobj`: pyfits object
                Science file, for which a detector to image correction 
                is available
                
        """
        assert isinstance(fobj, pyfits.NP_pyfits.HDUList)
        
        d2imfile = fobj[0].header['D2IMFILE']
        axiscorr = cls.getAxisCorr(d2imfile)
        d2imerr = pyfits.getdata(d2imfile, ext=1).max()
        if axiscorr == None:
            new_kw = {}
        else:
            new_kw = {'D2IMEXT': d2imfile, 'AXISCORR': axiscorr, 'D2IMERR': d2imerr}
            cls.applyDet2ImCorr(fobj, axiscorr)
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
        d2imfile = fileutil.osfn(fobj[0].header['D2IMFILE'])
        hdu = cls.createDgeoHDU(d2imfile, axiscorr)
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
    
    def createDgeoHDU(cls, d2imfile, axiscorr):
        
        d2im_data = pyfits.getdata(d2imfile, ext=1)
        d2im_hdr = cls.createDet2ImHdr(d2im_data.shape, axiscorr)
        hdu = pyfits.ImageHDU(header=d2im_hdr, data=d2im_data)
        
        return hdu
    
    createDgeoHDU = classmethod(createDgeoHDU)
    
    def createDet2ImHdr(cls, data_shape, axiscorr):
        """
        Creates a header for the D2IMARR extension based on the 
        reference file recorded in D2IMFILE keyword in the primary header.
        """
        
        naxis1 = data_shape[0]
        naxis2 = 0
        crpix1 = 0.0
        crpix2 = 0.0
        cdelt1 = 1.0
        cdelt2 = 1.0
        crval1 = 0.0
        crval2 = 0.0
        keys = ['XTENSION','BITPIX','NAXIS','NAXIS1','NAXIS2',
              'EXTNAME','EXTVER','PCOUNT','GCOUNT','CRPIX1',
                        'CDELT1','CRVAL1','CRPIX2','CDELT2','CRVAL2', 'AXISCORR']
                        
        comments = {'XTENSION': 'Image extension',
                    'BITPIX': 'IEEE floating point',
                    'NAXIS': 'Number of axes',
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
                    'CRVAL2': 'Image array pixel coordinate',
                    'AXISCORR': 'Direction in which the det2im correction is applied'}
        
        values = {'XTENSION': 'IMAGE',
                'BITPIX': -32,
                'NAXIS': 1,
                'NAXIS1': naxis1,
                'NAXIS2': naxis2,
                'EXTNAME': 'D2IMARR',
                'EXTVER':  1,
                'PCOUNT': 0,
                'GCOUNT': 1,
                'CRPIX1': crpix1,
                'CDELT1': cdelt1,
                'CRVAL1': crval1,
                'CRPIX2': crpix2,
                'CDELT2': cdelt2,
                'CRVAL2': crval2,
                'AXISCORR': axiscorr
                }
                    
        
        cdl = pyfits.CardList()
        for c in keys:
            cdl.append(pyfits.Card(key=c, value=values[c], comment=comments[c]))

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
    