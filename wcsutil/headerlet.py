from __future__ import division
import os
import tarfile
import tempfile
import time
import numpy as np
import warnings
import pyfits

from cStringIO import StringIO

from hstwcs import HSTWCS
import altwcs
from mappings import basic_wcs

import logging
logger = logging.getLogger(__name__)

def isWCSIdentical(scifile, file2, verbose=False):
    """
    Compares the WCS solution of 2 files.

    Parameters
    ----------
    scifile: file1
    file2: file2
    verbose: False or a python logging level
             (one of 'INFO', 'DEBUG' logging levels)
             (an integer representing a logging level)

    Notes
    -----
    These can be 2 science observations or 2 headerlets
    or a science observation and a headerlet. The two files
    have the same WCS solution if the following are the same:

    - rootname/destim
    - primary WCS
    - SIP coefficients
    - NPOL distortion
    - D2IM correction
    - Velocity aberation

    """
    if verbose:
        logger.setLevel(verbose)
    else:
        logger.setLevel(100)

    logger.info("Starting isWCSIdentical: %s" % time.asctime())

    result = True
    numsci1 = max(countExt(scifile, 'SCI'), countExt(scifile, 'SIPWCS'))
    numsci2 = max(countExt(file2, 'SCI'), countExt(file2, 'SIPWCS'))

    if numsci1 == 0 or numsci2 == 0 or numsci1 != numsci2:
        logger.info("Number of SCI and SIPWCS extensions do not match.")
        result = False

    if getRootname(scifile) != getRootname(file2):
        logger.info('Rootnames do not match.')
        result = False
    try:
        extname1 = pyfits.getval(scifile, 'EXTNAME', ext=('SCI', 1))
    except KeyError:
        extname1 = 'SIPWCS'
    try:
        extname2 = pyfits.getval(file2, 'EXTNAME', ext=('SCI', 1))
    except KeyError:
        extname2 = 'SIPWCS'

    for i in range(1, numsci1 + 1):
        w1 = HSTWCS(scifile, ext=(extname1, i))
        w2 = HSTWCS(file2, ext=(extname2, i))
        if not (w1.wcs.crval == w2.wcs.crval).all() or \
            not (w1.wcs.crpix == w2.wcs.crpix).all()  or \
            not  (w1.wcs.cd == w2.wcs.cd).all() or \
            not (np.array(w1.wcs.ctype) == np.array(w2.wcs.ctype)).all():
            logger.info('Primary WCSs do not match')
            result = False
        if w1.sip or w2.sip:
            if (w2.sip and not w1.sip) or (w1.sip and not w2.sip) or \
                not (w1.sip.a == w2.sip.a).all() or \
                not (w1.sip.b == w2.sip.b).all():
                logger.info('SIP coefficients do not match')
                result = False
        if w1.cpdis1 or w2.cipdis1:
            if w1.cpdis1 and not w2.cpdis1 or \
                w2.cpdis1 and not w1.cpdis1 or \
                not (w1.cpdis1.data == w2.cpdis1.data).all():
                logger.info('NPOL distortions do not match')
                result = False
        if w1.cpdis2 or w2.cpdis2:
            if w1.cpdis2 and not w2.cpdis2 or \
                w2.cpdis2 and not w1.cpdis2 or \
                not (w1.cpdis2.data == w2.cpdis2.data).all():
                logger.info('NPOL distortions do not match')
                result = False
        if w1.det2im1 or w2.det2im1:
            if w1.det2im1 and not w2.det2im1 or \
                w2.det2im1 and not w1.det2im1 or\
                not (w1.det2im1.data == w2.det2im1.data).all():
                logger.info('Det2Im corrections do not match')
                result =  False
        if w1.det2im2 or w2.det2im2:
            if w1.det2im2 and not w2.det2im2 or \
                w2.det2im2 and not w1.det2im2 or\
                not (w1.det2im2.data == w2.det2im2.data).all():
                logger.info('Det2Im corrections do not match')
                result = False
        if w1.vafactor != w2.vafactor:
            logger.info('VA factors do not match')
            result = False

    return result


# TODO: It would be logical for this to be part of the Headerlet class, perhaps
# as a classmethod
def createHeaderlet(fname, hdrname, destim=None, output=None, verbose=False):
    """
    Create a headerlet from a science observation

    Parameters
    ----------
    fname: string
           Name of file with science observation
    hdrname: string
           Name for the headerlet, stored in the primary header of the headerlet
    destim: string
           Destination image, stored in the primary header of the headerlet.
           If None ROOTNAME is used of the science observation is used.
           ROOTNAME has precedence, destim is used for observations without
           ROOTNAME in the primary header
    output: string
           Save the headerlet to the given filename.
    verbose: False or a python logging level
             (one of 'INFO', 'DEBUG' logging levels)
             (an integer representing a logging level)
    """

    if verbose:
        logger.setLevel(verbose)
    else:
        logger.setLevel(100)

    logger.info("Starting createHeaderlet: %s" % time.asctime())
    fmt="%Y-%m-%dT%H:%M:%S"
    phdukw = {'IDCTAB': True,
            'NPOLFILE': True,
            'D2IMFILE': True}
    fobj = pyfits.open(fname)
    if destim is None:
        try:
            destim = fobj[0].header['ROOTNAME']
        except KeyError:
            logger.exception('Required keyword "DESTIM" not found')
            print 'Please provide a value for the DESTIM keyword'
            raise
    if hdrname is None:
        logger.critical("Required keyword 'HDRNAME' not given")
        raise ValueError("Please provide a name for the headerlet, HDRNAME is "
                         "a required parameter.")



    altkeys = altwcs.wcskeys(fobj[('SCI', 1)].header)
    try:
        upwcsver = fobj[0].header.ascard['STWCSVER']
    except KeyError:
        upwcsver = pyfits.Card("STWCSVER", " ",
                               "Version of STWCS used to update the WCS")
    if 'O' in altkeys:
        altkeys.remove('O')
    numsci = countExt(fname, 'SCI')
    logger.debug("Number of 'SCI' extensions in file %s is %s"
                 % (fname, numsci))
    hdul = pyfits.HDUList()
    phdu = pyfits.PrimaryHDU()
    phdu.header.update('DESTIM', destim,
                       comment='Destination observation root name')
    phdu.header.update('HDRNAME', hdrname, comment='Headerlet name')
    phdu.header.update('DATE', time.strftime(fmt),
                       comment='Date FITS file was generated')
    phdu.header.ascard.append(upwcsver)

    updateRefFiles(fobj[0].header.ascard, phdu.header.ascard)
    phdu.header.update(key='VAFACTOR',
                       value=fobj[('SCI',1)].header.get('VAFACTOR', 1.))
    hdul.append(phdu)

    for e in range(1, numsci + 1):
        hwcs = HSTWCS(fname, ext=('SCI', e))
        h = hwcs.wcs2header(sip2hdr=True).ascard
        for ak in altkeys:
            awcs = HSTWCS(fname,ext=('SCI', e), wcskey=ak)
            h.extend(awcs.wcs2header(idc2hdr=False).ascard)
        h.insert(0, pyfits.Card(key='EXTNAME', value='SIPWCS',
                                comment='Extension name'))
        h.insert(1, pyfits.Card(key='EXTVER', value=e,
                                comment='Extension version'))

        fhdr = fobj[('SCI', e)].header.ascard
        if phdukw['NPOLFILE']:

            cpdis = fhdr['CPDIS*...']
            for c in range(1, len(cpdis) + 1):
                h.append(cpdis[c - 1])
                dp = fhdr['DP%s.*...' % c]
                h.extend(dp)

                try:
                    h.append(fhdr['CPERROR%s' % c])
                except KeyError:
                    pass

            try:
                h.append(fhdr['NPOLEXT'])
            except KeyError:
                pass

        if phdukw['D2IMFILE']:
            try:
                h.append(fhdr['D2IMEXT'])
            except KeyError:
                pass

            try:
                h.append(fhdr['AXISCORR'])
            except KeyError:
                logger.exception("Required keyword AXISCORR was not found in "
                                 "%s['SCI',%d]" % (fname, e))
                raise

            try:
                h.append(fhdr['D2IMERR'])
            except KeyError:
                h.append(pyfits.Card(key='DPERROR', value=0,
                                     comment='Maximum error of D2IMARR'))

        hdu = pyfits.ImageHDU(header=pyfits.Header(h))
        # temporary fix for pyfits ticket # 48
        hdu._extver = e
        hdul.append(hdu)
    numwdvarr = countExt(fname, 'WCSDVARR')
    numd2im = countExt(fname, 'D2IMARR')
    for w in range(1, numwdvarr + 1):
        hdu = fobj[('WCSDVARR', w)].copy()
        # temporary fix for pyfits ticket # 48
        hdu._extver = w
        hdul.append(hdu)
    for d in range(1, numd2im + 1):
        hdu = fobj[('D2IMARR', d)].copy()
        # temporary fix for pyfits ticket # 48
        hdu._extver = d
        hdul.append(hdu)
    if output is not None:
        # write the headerlet to a file
        if not output.endswith('_hdr.fits'):
            output = output + '_hdr.fits'
        hdul.writeto(output, clobber=True)
    fobj.close()
    return Headerlet(hdul)

def applyHeaderlet(hdrfile, destfile, destim=None, hdrname=None,
                   createheaderlet=True, verbose=False):
    """
    Apply headerlet 'hdrfile' to a science observation 'destfile'

    Parameters
    ----------
    hdrfile: string
             Headerlet file
    destfile: string
             File name of science observation wyhose WCS solution willbe updated
    destim: string or None (default)
            ROOTNAME of destfile (default)
            This string will be written as the DESTIM keyword in the headerlet,
            created from the old WCS solution
            Normally it should be None, unless the science file is missing the
            ROOTNAME kw, in which case the name of the file (stripped of .fits)
            should be specified.
    hdrname: string or None (default)
            will be the value of the HDRNAME keyword in the headerlet.
            It's not  required if createheaderlet is False
    createheaderlet: boolean
            True (default): before updating, create a headerlet with the
            WCS old solution.
    verbose: False or a python logging level
             (one of 'INFO', 'DEBUG' logging levels)
             (an integer representing a logging level)
    """

    logger.info("Starting applyHeaderlet: %s" % time.asctime())
    hlet = Headerlet(hdrfile)
    hlet.apply(destfile, destim=destim, hdrname=hdrname,
               createheaderlet=createheaderlet)

def updateRefFiles(source, dest):
    """
    Update the reference files name in the primary header of 'dest'
    using values from 'source'

    Parameters
    ----------
    source: pyfits.Header.ascardlist
    dest:   pyfits.Header.ascardlist
    """

    logger.info("Updating reference files")
    phdukw = {'IDCTAB': True,
            'NPOLFILE': True,
            'D2IMFILE': True}

    try:
        wind = dest.index_of('HISTORY')
    except KeyError:
        wind = len(dest)
    for key in phdukw.keys():
        try:
            value = source[key]
            dest.insert(wind, value)
        except KeyError:
            # TODO: I don't understand what the point of this is.  Is it meant
            # for logging purposes?  Right now it isn't used.
            phdukw[key] = False

def getRootname(fname):
    """
    returns the value of ROOTNAME or DESTIM
    """

    try:
        rootname = pyfits.getval(fname, 'ROOTNAME')
    except KeyError:
        rootname = pyfits.getval(fname, 'DESTIM')
    return rootname

def mapFitsExt2HDUListInd(fname, extname):
    """
    Map FITS extensions with 'EXTNAME' to HDUList indexes.
    """

    if isinstance(fname, basestring):
        f = pyfits.open(fname)
    else:
        f = fname
    d = {}
    for hdu in f:
        # TODO: Replace calls to header.has_key() with `in header` once
        # pyfits refactoring branch is in production use
        if hdu.header.has_key('EXTNAME') and hdu.header['EXTNAME'] == extname:
            extver = hdu.header['EXTVER']
            d[(extname, extver)] = f.index_of((extname, extver))
    f.close()
    return d

def countExt(fname, extname):
    """
    Return the number of extensions with EXTNAME in the file.
    """

    if isinstance(fname, basestring):
        f = pyfits.open(fname)
    else:
        f = fname
    numext = 0
    for hdu in f:
        if hdu.header.has_key('EXTNAME') and hdu.header['EXTNAME'] == extname:
            numext+=1
    logger.debug("file %s has %s extensions with extname=%s"
                 % (f.filename(), numext, extname))
    return numext


class Headerlet(pyfits.HDUList):
    """
    A Headerlet class
    Ref: http://stsdas.stsci.edu/stsci_python_sphinxdocs/stwcs/headerlet_def.html
    """

    def __init__(self, fobj, wcskeys=[], mode='copyonwrite', verbose=False):
        """
        Parameters
        ----------
        fobj:  string
                Name of headerlet file, file-like object, a list of HDU
                instances, or an HDUList instance
        wcskeys: python list
                a list of wcskeys to be included in the headerlet
                created from the old WCS solution before the
                science file is updated. If empty: all alternate (if any)
                WCSs are copied to the headerlet.
        mode: string, optional
                Mode with which to open the given file object
        """

        logger.info("Creating a Headerlet object from wcskeys %s" % wcskeys)
        self.wcskeys = wcskeys
        if not isinstance(fobj, list):
            fobj = pyfits.open(fobj, mode=mode)

        super(Headerlet, self).__init__(fobj)
        self.fname = self.filename()
        self.hdrname = self[0].header["HDRNAME"]
        self.stwcsver = self[0].header.get("STWCSVER", "")
        self.destim = self[0].header["DESTIM"]
        self.idctab = self[0].header.get("IDCTAB", "")
        self.npolfile = self[0].header.get("NPOLFILE", "")
        self.d2imfile = self[0].header.get("D2IMFILE", "")
        self.vafactor = self[1].header.get("VAFACTOR", 1) #None instead of 1?
        self.d2imerr = 0
        self.axiscorr = 1

    def apply(self, dest, destim=None, hdrname=None, createheaderlet=True):
        """
        Apply this headerlet to a file.
        """

        self.hverify()
        if self.verify_dest(dest):
            if createheaderlet:
                # TODO: Currently this does nothing...the created headerlet
                # isn't used--it should be given an appropriate name and
                # appended to the end of the file.
                orig_headerlet = createHeaderlet(dest, hdrname, destim)
            fobj = pyfits.open(dest, mode='update')
            self._delDestWCS(fobj)
            updateRefFiles(self[0].header.ascard, fobj[0].header.ascard)
            numsip = countExt(self, 'SIPWCS')
            for idx in range(1, numsip + 1):
                fhdr = fobj[('SCI', idx)].header.ascard
                siphdr = self[('SIPWCS', idx)].header.ascard
                # a minimal attempt to get the position of the WCS keywords group
                # in the header by looking for the PA_APER kw.
                # at least make sure the WCS kw are written befir the HISTORY kw
                # if everything fails, append the kw to the header
                try:
                    wind = fhdr.index_of('PA_APER')
                except KeyError:
                    try:
                        wind = fhdr.index_of('HISTORY')
                    except KeyError:
                        wind = len(fhdr)
                logger.debug("Inserting WCS keywords at index %s" % wind)
                # TODO: Drop .keys() when refactored pyfits comes into use
                for k in siphdr.keys():
                    if k not in ['XTENSION', 'BITPIX', 'NAXIS', 'PCOUNT',
                                 'GCOUNT','EXTNAME', 'EXTVER', 'ORIGIN',
                                 'INHERIT', 'DATE', 'IRAF-TLM']:
                        fhdr.insert(wind, siphdr[k])
                    else:
                        pass

            #! Always attach these extensions last. Otherwise their headers may
            # get updated with the other WCS kw.
            numwdvar = countExt(self, 'WCSDVARR')
            numd2im = countExt(self, 'D2IMARR')
            for idx in range(1, numwdvar + 1):
                fobj.append(self[('WCSDVARR', idx)].copy())
            for idx in range(1, numd2im + 1):
                fobj.append(self[('D2IMARR', idx)].copy())
            fobj.close()
        else:
            logger.critical("Observation %s cannot be updated with headerlet "
                            "%s" % (dest, self.hdrname))
            print "Observation %s cannot be updated with headerlet %s" \
                  % (dest, self.hdrname)


    def hverify(self):
        self.verify()
        assert(self[0].header.has_key('DESTIM') and
               self[0].header['DESTIM'].strip())
        assert(self[0].header.has_key('HDRNAME') and
               self[0].header['HDRNAME'].strip())
        assert(self[0].header.has_key('STWCSVER'))


    def verify_dest(self, dest):
        """
        verifies that the headerlet can be applied to the observation

        DESTIM in the primary header of the headerlet must match ROOTNAME
        of the science file (or the name of the destination file)
        """

        try:
            droot = pyfits.getval(dest, 'ROOTNAME')
        except KeyError:
            logger.debug("Keyword 'ROOTNAME' not found in destination file")
            droot = dest.split('.fits')[0]
        if droot == self.destim:
            logger.debug("verify_destim() returned True")
            return True
        else:
            logger.debug("verify_destim() returned False")
            return False

    def tofile(self, fname, destim=None, hdrname=None, clobber=False):
        if not destim or not hdrname:
            self.hverify()
        self.writeto(fname, clobber=clobber)

    def _delDestWCS(self, dest):
        """
        Delete the WCS of a science file
        """

        logger.info("Deleting all WCSs of file %s"
                    % dest.fileinfo()['filename'])
        numext = len(dest)

        for idx in range(numext):
            self._removeD2IM(dest[idx])
            self._removeSIP(dest[idx])
            self._removeLUT(dest[idx])
            self._removePrimaryWCS(dest[idx])
            self._removeIDCCoeffs(dest[idx])
            try:
                del dest[idx].header.ascard['VAFACTOR']
            except KeyError:
                pass

        self._removeAltWCS(dest, ext=range(numext))
        numwdvarr = countExt(dest, 'WCSDVARR')
        numd2im = countExt(dest, 'D2IMARR')
        for idx in range(1, numwdvarr + 1):
            del dest[('WCSDVARR', idx)]
        for idx in range(1, numd2im + 1):
            del dest[('D2IMARR', idx)]

    def _removeSIP(self, ext):
        """
        Remove the SIP distortion of a FITS extension
        """

        logger.debug("Removing SIP distortion from (%s, %s)"
                     % (ext.name, ext._extver))
        for prefix in ['A', 'B', 'AP', 'BP']:
            try:
                order = ext.header[prefix + '_ORDER']
                del ext.header[prefix + '_ORDER']
            except KeyError:
                continue
            for i in range(order + 1):
                for j in range(order + 1):
                    key = prefix + '_%d_%d' % (i, j)
                    try:
                        del ext.header[key]
                    except KeyError:
                        pass
        try:
            del ext.header['IDCTAB']
        except KeyError:
            pass

    def _removeLUT(self, ext):
        """
        Remove the Lookup Table distortion of a FITS extension
        """

        logger.debug("Removing LUT distortion from (%s, %s)"
                     % (ext.name, ext._extver))
        try:
            cpdis = ext.header['CPDIS*']
        except KeyError:
            return
        try:
            for c in range(1, len(cpdis) + 1):
                del ext.header['DP%s.*...' % c]
                del ext.header[cpdis[c - 1].key]
            del ext.header['CPERR*']
            del ext.header['NPOLFILE']
            del ext.header['NPOLEXT']
        except KeyError:
            pass

    def _removeD2IM(self, ext):
        """
        Remove the Detector to Image correction of a FITS extension
        """

        logger.debug("Removing D2IM correction from (%s, %s)"
                     % (ext.name, ext._extver))
        d2imkeys = ['D2IMFILE', 'AXISCORR', 'D2IMEXT', 'D2IMERR']
        for k in d2imkeys:
            try:
                del ext.header[k]
            except KeyError:
                pass

    def _removeAltWCS(self, dest, ext):
        """
        Remove Alternate WCSs of a FITS extension.
        A WCS with wcskey 'O' is never deleted.
        """

        dkeys = altwcs.wcskeys(dest[('SCI', 1)].header)
        logger.debug("Removing alternate WCSs with keys %s from %s"
                     % (dkeys, dest.filename()))
        for k in dkeys:
            altwcs.deleteWCS(dest, ext=ext, wcskey=k)

    def _removePrimaryWCS(self, ext):
        """
        Remove the primary WCS of a FITS extension
        """

        logger.debug("Removing Primary WCS from (%s, %s)"
                     % (ext.name, ext._extver))
        naxis = ext.header.ascard['NAXIS'].value
        for key in basic_wcs:
            for i in range(1, naxis + 1):
                try:
                    del ext.header.ascard[key + str(i)]
                except KeyError:
                    pass
        try:
            del ext.header.ascard['WCSAXES']
        except KeyError:
            pass

    def _removeIDCCoeffs(self, ext):
        """
        Remove IDC coefficients of a FITS extension
        """

        logger.debug("Removing IDC coefficient from (%s, %s)"
                     % (ext.name, ext._extver))
        coeffs = ['OCX10', 'OCX11', 'OCY10', 'OCY11', 'IDCSCALE']
        for k in coeffs:
            try:
                del ext.header.ascard[k]
            except KeyError:
                pass


class HeaderletHDU(pyfits.core._NonstandardExtHDU):
    """
    A non-standard extension HDU for encapsulating Headerlets in a file.  These
    HDUs have an extension type of HDRLET and their EXTNAME is derived from the
    Headerlet's HDRNAME.

    The data itself is a tar file containing a single file, which is the
    Headerlet file itself.  The file name is derived from the HDRNAME keyword,
    and should be in the form `<HDRNAME>_hdr.fits`.  If the COMPRESS keyword
    evaluates to `True`, the tar file is compressed with gzip compression.

    The Headerlet contained in the HDU's data can be accessed by the
    `headerlet` attribute.
    """

    _xtn = _extension = 'HDRLET'

    def __init__(self, data=None, header=None):
        super(HeaderletHDU, self).__init__(data=data, header=header)
        # TODO: This can be removed after the next pyfits release, but for now
        # the _ExtensionHDU base class sets self._xtn = '' in its __init__().
        self._xtn = self._extension
        # For some reason _NonstandardExtHDU.__init__ sets self.name = None,
        # even if it's already been set by the EXTNAME keyword in
        # _ExtensionHDU.__init__() -_-;
        if header and header.has_key('EXTNAME') and not self.name:
            self.name = header['EXTNAME']
        # self._extver, if set, is still preserved.  From
        # _ExtensionHDU.__init__()...totally inconsistent.

    def __getattr__(self, attr):
        if attr == 'data':
            size = self.size()
            self._file.seek(self._datLoc)
            self.__dict__[attr] = self._file.read(size)
        elif attr == 'headerlet':
            self._file.seek(self._datLoc)
            s = StringIO()
            # Read the data into a StringIO--reading directly from the file
            # won't work (at least for gzipped files) due to problems deep
            # within the gzip module that make it difficult to read gzip files
            # embedded in another file
            s.write(self._file.read(self.size()))
            s.seek(0)
            if self._header['COMPRESS']:
                mode = 'r:gz'
            else:
                mode = 'r'
            t = tarfile.open(mode=mode, fileobj=s)
            members = t.getmembers()
            if not len(members):
                raise ValueError('The Headerlet contents are missing.')
            elif len(members) > 1:
                warnings.warn('More than one file is contained in this '
                              'only the headerlet file should be present.')
            hlt_name = self.name + '_hdr.fits'
            try:
                hlt_info = t.getmember(hlt_name)
            except KeyError:
                warnings.warn('The file %s was missing from the HDU data.  '
                              'Assuming that the first file in the data is '
                              'headerlet file.')
                hlt_info = members[0]
            hlt_file = t.extractfile(hlt_info)
            # hlt_file is a file-like object
            hlt = Headerlet(hlt_file, mode='readonly')
            self.__dict__[attr] = hlt
        else:
            return pyfits.core._ValidHDU.__getattr__(self, attr)
        try:
            return self.__dict__[attr]
        except KeyError:
            raise AttributeError(attr)

    @classmethod
    def fromheaderlet(cls, headerlet, compress=False):
        """
        Creates a new HeaderletHDU from a given Headerlet object.

        Parameters
        ----------
        headerlet : Headerlet
            A valid Headerlet object.
        compress : bool, optional
            Gzip compress the headerlet data.
        """

        phdu = headerlet[0]
        phduhdr = phdu.header
        hlt_filename = phdu.header['HDRNAME'] + '_hdr.fits'

        # TODO: As it stands there's no good way to write out an HDUList in
        # memory, since it automatically closes the given file-like object when
        # it's done writing.  I'd argue that if passed an open file handler it
        # should not close it, but for now we'll have to write to a temp file.
        fd, name = tempfile.mkstemp()
        try:
            f = os.fdopen(fd, 'rb+')
            headerlet.writeto(f)
            # The tar file itself we'll write in memory, as it should be
            # relatively small
            if compress:
                mode = 'w:gz'
            else:
                mode = 'w'
            s = StringIO()
            t = tarfile.open(mode=mode, fileobj=s)
            t.add(name, arcname=hlt_filename)
            t.close()
        finally:
            os.remove(name)

        cards = [
            pyfits.Card('XTENSION', cls._extension, 'Headerlet extension'),
            pyfits.Card('BITPIX',    8, 'array data type'),
            pyfits.Card('NAXIS',     1, 'number of array dimensions'),
            pyfits.Card('NAXIS1',    len(s.getvalue()), 'Axis length'),
            pyfits.Card('PCOUNT',    0, 'number of parameters'),
            pyfits.Card('GCOUNT',    1, 'number of groups'),
            pyfits.Card('EXTNAME', phdu.header['HDRNAME'],
                        'name of the headerlet extension'),
            phdu.header.ascard['HDRNAME'],
            phdu.header.ascard['DATE'],
            pyfits.Card('SIPNAME', headerlet['SIPWCS', 1].header['WCSNAMEA'],
                        'SIP distortion model name'),
            phdu.header.ascard['NPOLFILE'],
            phdu.header.ascard['D2IMFILE'],
            pyfits.Card('COMPRESS', compress, 'Uses gzip compression')
        ]

        header = pyfits.Header(pyfits.CardList(cards))

        hdu = cls(data=pyfits.DELAYED, header=header)
        hdu._file = s
        hdu._datLoc = 0
        return hdu

    @classmethod
    def match_header(cls, header):
        """
        This is a class method used in the pyfits refactoring branch to
        recognize whether or not this class should be used for instantiating
        an HDU object based on values in the header.

        It is included here for forward-compatibility.
        """

        card = header.ascard[0]
        if card.key != 'XTENSION':
            return False
        xtension = card.value.rstrip()
        return xtension == cls._extension

    # TODO: Add header verification

    def _summary(self):
        # TODO: Perhaps make this more descriptive...
        return '%-10s  %-12s  %4d' % (self.name, self.__class__.__name__,
                                      len(self._header.ascard))

# Monkey-patch pyfits to add minimal support for HeaderletHDUs
# TODO: Get rid of this ASAP!!! (it won't be necessary with the pyfits
# refactoring branch)
__old_updateHDUtype = pyfits.Header._updateHDUtype
def __updateHDUtype(self):
    if HeaderletHDU.match_header(self):
        self._hdutype = HeaderletHDU
    else:
        __old_updateHDUtype(self)
pyfits.Header._updateHDUtype = __updateHDUtype