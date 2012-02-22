from __future__ import division # confidence high

import time
import pyfits
from stsci.tools import fileutil
import utils

import logging
logger = logging.getLogger('stwcs.updatewcs.Det2IM')

class DET2IMCorr(object):
    """
    Stores a small correction to the detector coordinates as a d2imarr
    extension in the science file.

    Notes
    -----
    For the case of ACS/WFC every 68th column is wider than the rest.
    To compensate for this a small correction needs to be applied to the
    detector coordinates. We call this a detector to image transformation.
    The so obtained image coordinates are the input to all other distortion
    corrections. The correction is originally stored in an external
    reference file pointed to by 'D2IMFILE' keyword in the primary header.
    This class attaches the correction array as an extension to the science
    file with extname = `d2imarr`.

    Other keywords used in this correction are:

    `AXISCORR`: integer (1 or 2) - axis to which the detector to image
                correction is applied

    `D2IMEXT`:  string = name of reference file which was used to create
                the lookup table extension

    `D2IMERR`:  float, optional - maximum value of the correction

    """
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
        logger.info("\n\tStarting Det2IM Correction: %s" % time.asctime())
        try:
            assert isinstance(fobj, pyfits.HDUList)
        except AssertionError:
            logger.exception('\n\tInput must be a pyfits.HDUList object')
            raise

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
            direction = pyfits.getval(refname, 'EXTNAME', ext=1)
            if direction == 'DX': return 1
            elif direction == 'DY': return 2
            else:
                logger.warning('\n\tDET2IM correction expects the reference file to have \
                an EXTNAME keyword of value "DX"  or "DY", EXTNAMe %s detected' % direction)
                return None
        except KeyError:
            logger.exception("\n\tD2IMFILE %s is missing EXTNAME keyword. Unable to determine axis \
            to which to apply the correction." % refname)
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
        d2im = pyfits.open(d2imfile)
        data_shape = d2im[1].shape
        naxis = d2im[1].header['NAXIS']
        d2im_phdr = d2im[0].header
        d2im.close()

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


        cdl = []
        for key in kw_comm0.keys():
            cdl.append((key, kw_val0[key], kw_comm0[key]))
        for key in kw_comm1.keys():
            cdl.append((key, kw_val1[key], kw_comm1[key]))
        # Now add keywords from NPOLFILE header to document source of calibration
        # include all keywords after and including 'FILENAME' from header
        start_indx = -1
        end_indx = 0
        for i, c in enumerate(d2im_phdr):
            if c == 'FILENAME':
                start_indx = i
            if c == '': # remove blanks from end of header
                end_indx = i+1
                break
        if start_indx >= 0:
            for card in d2im_phdr.cards[start_indx:end_indx]:
                cdl.append(card)

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

