import logging
import time

import numpy as np
from astropy.io import fits

from stsci.tools import fileutil

logger = logging.getLogger('stwcs.updatewcs.npol')


class NPOLCorr(object):
    """
    Defines a Lookup table prior distortion correction as per WCS paper IV.
    It uses a reference file defined by the NPOLFILE (suffix 'NPL') keyword
    in the primary header.

    Notes
    -----
    - Using extensions in the reference file create a WCSDVARR extensions
      and add them to the science file.
    - Add record-valued keywords to the science extension header to describe
      the lookup tables.
    - Add a keyword 'NPOLEXT' to the science extension header to store
      the name of the reference file used to create the WCSDVARR extensions.

    If WCSDVARR extensions exist and `NPOLFILE` is different from `NPOLEXT`,
    a subsequent update will overwrite the existing extensions.
    If WCSDVARR extensions were not found in the science file, they will be added.

    It is assumed that the NPL reference files were created to work with IDC tables
    but will be applied with SIP coefficients. A transformation is applied to correct
    for the fact that the lookup tables will be applied before the first order coefficients
    which are in the CD matrix when the SIP convention is used.
    """

    def updateWCS(cls, fobj):
        """
        Parameters
        ----------
        fobj : `astropy.io.fits.HDUList` object
            Science file, for which a distortion correction in a NPOLFILE is available

        """
        logger.info("\n\tStarting NPOL: %s" % time.asctime())
        try:
            assert isinstance(fobj, fits.HDUList)
        except AssertionError:
            logger.exception('\n\tInput must be a fits.HDUList object')
            raise

        cls.applyNPOLCorr(fobj)
        nplfile = fobj[0].header['NPOLFILE']

        new_kw = {'NPOLEXT': nplfile}
        return new_kw

    updateWCS = classmethod(updateWCS)

    def applyNPOLCorr(cls, fobj):
        """
        For each science extension in a fits file object:
            - create a WCSDVARR extension
            - update science header
            - add/update NPOLEXT keyword
        """
        nplfile = fileutil.osfn(fobj[0].header['NPOLFILE'])
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
                header = ext.header
                # get the data arrays from the reference file and transform
                # them for use with SIP
                dx, dy = cls.getData(nplfile, ccdchip)
                idccoeffs = cls.getIDCCoeffs(header)

                if idccoeffs is not None:
                    dx, dy = cls.transformData(dx, dy, idccoeffs)

                # Determine EXTVER for the WCSDVARR extension from the
                # NPL file (EXTNAME, EXTVER) kw.
                # This is used to populate DPj.EXTVER kw
                wcsdvarr_x_version = 2 * extversion - 1
                wcsdvarr_y_version = 2 * extversion
                for ename in zip(['DX', 'DY'], [wcsdvarr_x_version, wcsdvarr_y_version], [dx, dy]):
                    error_val = ename[2].max()
                    cls.addSciExtKw(header, wdvarr_ver=ename[1], npol_extname=ename[0], error_val=error_val)
                    hdu = cls.createNpolHDU(header, npolfile=nplfile,
                                            wdvarr_ver=ename[1], npl_extname=ename[0],
                                            data=ename[2], ccdchip=ccdchip)
                    if wcsdvarr_ind:
                        fobj[wcsdvarr_ind[ename[1]]] = hdu
                    else:
                        fobj.append(hdu)

    applyNPOLCorr = classmethod(applyNPOLCorr)

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
        logger.debug("A map of WSCDVARR externsions %s" % wcsd)
        return wcsd

    getWCSIndex = classmethod(getWCSIndex)

    def addSciExtKw(cls, hdr, wdvarr_ver=None, npol_extname=None, error_val=0.0):
        """
        Adds kw to sci extension to define WCSDVARR lookup table extensions

        """
        if npol_extname == 'DX':
            j = 1
        else:
            j = 2

        cperror = 'CPERR%s' % j
        cpdis = 'CPDIS%s' % j
        dpext = 'DP%s.' % j + 'EXTVER'
        dpnaxes = 'DP%s.' % j + 'NAXES'
        dpaxis1 = 'DP%s.' % j + 'AXIS.1'
        dpaxis2 = 'DP%s.' % j + 'AXIS.2'
        keys = [cperror, cpdis, dpext, dpnaxes, dpaxis1, dpaxis2]
        values = {cperror: error_val,
                  cpdis: 'Lookup',
                  dpext: wdvarr_ver,
                  dpnaxes: 2,
                  dpaxis1: 1,
                  dpaxis2: 2}

        comments = {cperror: 'Maximum error of NPOL correction for axis %s' % j,
                    cpdis: 'Prior distortion function type',
                    dpext: 'Version number of WCSDVARR extension containing lookup distortion table',
                    dpnaxes: 'Number of independent variables in distortion function',
                    dpaxis1: 'Axis number of the jth independent variable in a distortion function',
                    dpaxis2: 'Axis number of the jth independent variable in a distortion function'
                    }
        # Look for HISTORY keywords. If present, insert new keywords before them
        before_key = 'HISTORY'
        if before_key not in hdr:
            before_key = None

        for key in keys:
            hdr.set(key, value=values[key], comment=comments[key], before=before_key)

    addSciExtKw = classmethod(addSciExtKw)

    def getData(cls, nplfile, ccdchip):
        """
        Get the data arrays from the reference NPOL files
        Make sure 'CCDCHIP' in the npolfile matches "CCDCHIP' in the science file.
        """
        npl = fits.open(nplfile)
        for ext in npl:
            nplextname  = ext.header.get('EXTNAME', "")
            nplccdchip  = ext.header.get('CCDCHIP', 1)
            if nplextname == 'DX' and nplccdchip == ccdchip:
                xdata = ext.data.copy()
                continue
            elif nplextname == 'DY' and nplccdchip == ccdchip:
                ydata = ext.data.copy()
                continue
            else:
                continue
        npl.close()
        return xdata, ydata
    getData = classmethod(getData)

    def transformData(cls, dx, dy, coeffs):
        """
        Transform the NPOL data arrays for use with SIP
        """
        ndx, ndy = np.dot(coeffs, [dx.ravel(), dy.ravel()]).astype(np.float32)
        ndx.shape = dx.shape
        ndy.shape = dy.shape
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
            coeffs = np.array([[ocx11, ocx10], [ocy11, ocy10]], dtype=np.float32)
        except KeyError:
            logger.exception('\n\tFirst order IDCTAB coefficients are not available. \n\
            Cannot convert SIP to IDC coefficients.')
            return None
        try:
            idcscale = header['IDCSCALE']
        except KeyError:
            logger.exception("IDCSCALE not found in header - setting it to 1.")
            idcscale = 1

        return np.linalg.inv(coeffs / idcscale)

    getIDCCoeffs = classmethod(getIDCCoeffs)

    def createNpolHDU(cls, sciheader, npolfile=None, wdvarr_ver=1, npl_extname=None,
                      data=None, ccdchip=1):
        """
        Creates an HDU to be added to the file object.
        """
        hdr = cls.createNpolHdr(sciheader, npolfile=npolfile, wdvarr_ver=wdvarr_ver,
                                npl_extname=npl_extname, ccdchip=ccdchip)
        hdu = fits.ImageHDU(header=hdr, data=data)
        return hdu

    createNpolHDU = classmethod(createNpolHDU)

    def createNpolHdr(cls, sciheader, npolfile, wdvarr_ver, npl_extname, ccdchip):
        """
        Creates a header for the WCSDVARR extension based on the NPOL reference file
        and sci extension header. The goal is to always work in image coordinates
        (also for subarrays and binned images. The WCS for the WCSDVARR extension
        i ssuch that a full size npol table is created and then shifted or scaled
        if the science image is a subarray or binned image.
        """
        npl = fits.open(npolfile)
        npol_phdr = npl[0].header
        for ext in npl:
            try:
                nplextname = ext.header['EXTNAME']
                nplextver = ext.header['EXTVER']
            except KeyError:
                continue
            nplccdchip = cls.get_ccdchip(npl, extname=nplextname, extver=nplextver)
            if nplextname == npl_extname and nplccdchip == ccdchip:
                npol_header = ext.header
                break
            else:
                continue
        npl.close()

        naxis = npl[1].header['NAXIS']
        ccdchip = nplextname  # npol_header['CCDCHIP']

        kw = {'NAXIS': 'Size of the axis',
              'CDELT': 'Coordinate increment along axis',
              'CRPIX': 'Coordinate system reference pixel',
              'CRVAL': 'Coordinate system value at reference pixel',
              }

        kw_comm1 = {}
        kw_val1 = {}
        for key in kw.keys():
            for i in range(1, naxis + 1):
                si = str(i)
                kw_comm1[key + si] = kw[key]

        for i in range(1, naxis + 1):
            si = str(i)
            kw_val1['NAXIS' + si] = npol_header.get('NAXIS' + si)
            kw_val1['CDELT' + si] = npol_header.get('CDELT' + si, 1.0) * \
                sciheader.get('LTM' + si + '_' + si, 1)
            kw_val1['CRPIX' + si] = npol_header.get('CRPIX' + si, 0.0)
            kw_val1['CRVAL' + si] = (npol_header.get('CRVAL' + si, 0.0) -
                                     sciheader.get('LTV' + str(i), 0))

        kw_comm0 = {'XTENSION': 'Image extension',
                    'BITPIX': 'IEEE floating point',
                    'NAXIS': 'Number of axes',
                    'EXTNAME': 'WCS distortion array',
                    'EXTVER': 'Distortion array version number',
                    'PCOUNT': 'Special data area of size 0',
                    'GCOUNT': 'One data group',
                    }

        kw_val0 = {'XTENSION': 'IMAGE',
                   'BITPIX': -32,
                   'NAXIS': naxis,
                   'EXTNAME': 'WCSDVARR',
                   'EXTVER': wdvarr_ver,
                   'PCOUNT': 0,
                   'GCOUNT': 1,
                   'CCDCHIP': ccdchip,
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
        for i, c in enumerate(npol_phdr):
            if c == 'FILENAME':
                start_indx = i
            if c == '':  # remove blanks from end of header
                end_indx = i + 1
                break
        if start_indx >= 0:
            for card in npol_phdr.cards[start_indx: end_indx]:
                cdl.append(card)

        hdr = fits.Header(cards=cdl)

        return hdr

    createNpolHdr = classmethod(createNpolHdr)

    def get_ccdchip(cls, fobj, extname, extver):
        """
        Given a science file or npol file determine CCDCHIP
        """
        ccdchip = 1
        if fobj[0].header['INSTRUME'] == 'ACS' and fobj[0].header['DETECTOR'] == 'WFC':
            ccdchip = fobj[extname, extver].header['CCDCHIP']
        elif fobj[0].header['INSTRUME'] == 'WFC3' and fobj[0].header['DETECTOR'] == 'UVIS':
            ccdchip = fobj[extname, extver].header['CCDCHIP']
        elif fobj[0].header['INSTRUME'] == 'WFPC2':
            ccdchip = fobj[extname, extver].header['DETECTOR']
        elif fobj[0].header['INSTRUME'] == 'STIS':
            ccdchip = fobj[extname, extver].header['DETECTOR']
        elif fobj[0].header['INSTRUME'] == 'NICMOS':
            ccdchip = fobj[extname, extver].header['CAMERA']
        return ccdchip

    get_ccdchip = classmethod(get_ccdchip)
