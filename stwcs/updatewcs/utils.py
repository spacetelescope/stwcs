import os
from astropy.io import fits
from stsci.tools import fileutil

import logging
logger = logging.getLogger("stwcs.updatewcs.utils")

def get_file(input):
    """Return FITSObject object for input file which matches user input"""
    if isinstance(input, FITSObject):
        return input

    if isinstance(input, str):
        obj = FileObject(input)
    elif isinstance(input, fits.HDUList):
        obj = HDUObject(input)
    else:
        msg = "Unrecognized input type: {}".format(type(input))
        raise ValueError(msg)
    return obj

class FITSObject:
    """Abstract class to isolate file operations from files or in-memory objects"""
    def __init__(self, input):
        pass

    def open(self):
        """Return HDUList object"""
        pass

    def getval(self, kw, **kwargs):
        """Return keyword value for specified keyword"""
        pass

    def getheader(self):
        """Return requested header from file/object"""
        pass

class FileObject(FITSObject):
    def __init__(self, input):
        FITSObject.__init__(self, input)
        self.name = input
        self.close_hdu = True

    def open(self, mode='readonly', memmap=None, save_backup=False, cache=True, lazy_load_hdus=None, **kwargs):
        """Pass-through interface for fits.open"""
        self.hdulist = fits.open(self.name, mode=mode, memmap=memmap,
                            save_backup=save_backup, cache=cache,
                            lazy_load_hdus=lazy_load_hdus, **kwargs)
        return self.hdulist

    def close(self):
        self.hdulist.close()

    def getval(self,kw, *args, **kwargs):
        val = fits.getval(self.name, kw, *args, **kwargs)
        return val

    def getheader(self, *args, **kwargs):
        hdr = fits.getheader(self.name, *args, **kwargs)
        return hdr

class HDUObject(FITSObject):
    def __init__(self, input):
        FITSObject.__init__(self, input)
        try:
            assert isinstance(input, fits.HDUList)
        except AssertionError:
            msg = "HDUObject can only be initialized using a FITS.HDUList object"
            raise ValueError(msg)
        self.hdulist = input

        # Set 'filename' for object
        if 'rootname' in self.hdulist[0].header:
            self.name = self.hdulist[0].header['rootname']
        else:
            self.name = None

        self.close_hdu = False


    def open(self, *args, **kwargs):
        mode = kwargs.get('mode', 'update')
        try:
            assert (mode == self.hdulist._file.mode)
        except AssertionError:
            msg = 'Input HDUList not opened in "update" mode!'
            raise ValueError(msg)
        return self.hdulist

    def close(self):
        pass

    def getval(self, kw, *args, **kwargs):
        extn = kwargs.get('ext', 0)
        val = self.hdulist[extn].header[kw]
        return val

    def getheader(self, *args, **kwargs):
        extn = kwargs.get('ext', 0)
        hdr = self.hdulist[extn].header
        return hdr


def diff_angles(a, b):
    """
    Perform angle subtraction a-b taking into account
    small-angle differences across 360degree line.
    """

    diff = a - b

    if diff > 180.0:
        diff -= 360.0

    if diff < -180.0:
        diff += 360.0

    return diff


def getBinning(fobj, extver=1):
    # Return the binning factor
    binned = 1
    if fobj[0].header['INSTRUME'] == 'WFPC2':
        mode = fobj[0].header.get('MODE', "")
        if mode == 'AREA': binned = 2
    else:
        binned = fobj['SCI', extver].header.get('BINAXIS', 1)
    return binned


def updateNEXTENDKw(fobj):
    """
    Updates PRIMARY header with correct value for NEXTEND, if present.

    Parameters
    -----------
    fobj : `astropy.io.fits.HDUList`
        The FITS object for file opened in `update` mode

    """
    if 'nextend' in fobj[0].header:
        fobj[0].header['nextend'] = len(fobj) - 1


def extract_rootname(kwvalue, suffix=""):
    """ Returns the rootname from a full reference filename

        If a non-valid value (any of ['','N/A','NONE','INDEF',None]) is input,
            simply return a string value of 'NONE'

        This function will also replace any 'suffix' specified with a blank.
    """
    # check to see whether a valid kwvalue has been provided as input
    if kwvalue.strip() in ['', 'N/A', 'NONE', 'INDEF', None]:
        return 'NONE'  # no valid value, so return 'NONE'

    # for a valid kwvalue, parse out the rootname
    # strip off any environment variable from input filename, if any are given
    if '$' in kwvalue:
        fullval = kwvalue[kwvalue.find('$') + 1: ]
    else:
        fullval = kwvalue
    # Extract filename without path from kwvalue
    fname = os.path.basename(fullval).strip()

    # Now, rip out just the rootname from the full filename
    rootname = fileutil.buildNewRootname(fname)

    # Now, remove any known suffix from rootname
    rootname = rootname.replace(suffix, '')
    return rootname.strip()


def construct_distname(fobj, wcsobj):
    """
    This function constructs the value for the keyword 'DISTNAME'.
    It relies on the reference files specified by the keywords 'IDCTAB',
    'NPOLFILE', and 'D2IMFILE'.

    The final constructed value will be of the form:
        <idctab rootname>-<npolfile rootname>-<d2imfile rootname>
    and have a value of 'NONE' if no reference files are specified.
    """
    idcname = extract_rootname(fobj[0].header.get('IDCTAB', "NONE"),
                               suffix='_idc')
    if (idcname is None or idcname == 'NONE') and wcsobj.sip is not None:
        idcname = 'UNKNOWN'

    npolname, npolfile = build_npolname(fobj)
    if npolname is None and wcsobj.cpdis1 is not None:
        npolname = 'UNKNOWN'

    d2imname, d2imfile = build_d2imname(fobj)
    if d2imname is None and wcsobj.det2im is not None:
        d2imname = 'UNKNOWN'

    sipname, idctab = build_sipname(fobj)

    distname = build_distname(sipname, npolname, d2imname)
    return {'DISTNAME': distname, 'SIPNAME': sipname}


def build_distname(sipname, npolname, d2imname):
    """
    Core function to build DISTNAME keyword value without the HSTWCS input.
    """
    distname = sipname.strip()
    if npolname != 'NONE' or d2imname != 'NONE':
        if d2imname != 'NONE':
            distname += '-' + npolname.strip() + '-' + d2imname.strip()
        else:
            distname += '-' + npolname.strip()

    return distname


def build_default_wcsname(idctab):
    idcname = extract_rootname(idctab, suffix='_idc')
    wcsname = 'IDC_' + idcname
    return wcsname


def build_sipname(fobj, fname=None, sipname=None):
    """
    Build a SIPNAME from IDCTAB

    Parameters
    ----------
    fobj : `astropy.io.fits.HDUList`
        file object
    fname : string
          science file name (to be used if ROOTNAME is not present
    sipname : string
          user supplied SIPNAME keyword

    Returns
    -------
    sipname, idctab
    """
    try:
        idctab = fobj[0].header['IDCTAB']
        idcname = extract_rootname(idctab, suffix='_idc')
    except KeyError:
        idctab = 'N/A'
        idcname = 'N/A'
    if not fname:
        try:
            fname = fobj.filename()
        except:
            fname = " "
    if not sipname:
        try:
            sipname = fobj[0].header["SIPNAME"]
            if idcname == 'N/A' or idcname not in sipname:
                raise KeyError
        except KeyError:
            if idcname != 'N/A':
                try:
                    rootname = fobj[0].header['rootname']
                except KeyError:
                    rootname = fname
                sipname = rootname + '_' + idcname
            else:
                if 'A_ORDER' in fobj[1].header or 'B_ORDER' in fobj[1].header:
                    sipname = 'UNKNOWN'
                else:
                    idcname = 'NOMODEL'

    return sipname, idctab


def build_npolname(fobj, npolfile=None):
    """
    Build a NPOLNAME from NPOLFILE

    Parameters
    ----------
    fobj : `astropy.io.fits.HDUList`
        file object
    npolfile : string
        user supplied NPOLFILE keyword

    Returns
    -------
    npolname, npolfile
    """
    if not npolfile:
        try:
            npolfile = fobj[0].header["NPOLFILE"]
        except KeyError:
            npolfile = ' '
            if fileutil.countExtn(fobj, 'WCSDVARR'):
                npolname = 'UNKNOWN'
            else:
                npolname = 'NOMODEL'
        npolname = extract_rootname(npolfile, suffix='_npl')
        if npolname == 'NONE':
            npolname = 'NOMODEL'
    else:
        npolname = extract_rootname(npolfile, suffix='_npl')
        if npolname == 'NONE':
            npolname = 'NOMODEL'
    return npolname, npolfile


def build_d2imname(fobj, d2imfile=None):
    """
    Build a D2IMNAME from D2IMFILE

    Parameters
    ----------
    fobj : `astropy.io.fits.HDUList`
        file object
    d2imfile : string
        user supplied NPOLFILE keyword

    Returns
    -------
    d2imname, d2imfile
    """
    if not d2imfile:
        try:
            d2imfile = fobj[0].header["D2IMFILE"]
        except KeyError:
            d2imfile = 'N/A'
            if fileutil.countExtn(fobj, 'D2IMARR'):
                d2imname = 'UNKNOWN'
            else:
                d2imname = 'NOMODEL'
        d2imname = extract_rootname(d2imfile, suffix='_d2i')
        if d2imname == 'NONE':
            d2imname = 'NOMODEL'
    else:
        d2imname = extract_rootname(d2imfile, suffix='_d2i')
        if d2imname == 'NONE':
            d2imname = 'NOMODEL'
    return d2imname, d2imfile


def remove_distortion(fname, dist_keyword):
    logger.info("Removing distortion {0} from file {0}".format(dist_keyword, fname))
    from ..wcsutil import altwcs
    if dist_keyword == "NPOLFILE":
        extname = "WCSDVARR"
        keywords = ["CPERR*", "DP1.*", "DP2.*", "CPDIS*", "NPOLEXT"]
    elif dist_keyword == "D2IMFILE":
        extname = "D2IMARR"
        keywords = ["D2IMERR*", "D2IM1.*", "D2IM2.*", "D2IMDIS*", "D2IMEXT"]
    else:
        raise AttributeError("Unrecognized distortion keyword "
                             "{0} when attempting to remove distortion".format(dist_keyword))
    if isinstance(fname, FITSObject):
        f = fname.open(mode='update')
    else:
        f = fits.open(fname, mode="update")
    ext_mapping = altwcs.mapFitsExt2HDUListInd(f, "SCI").values()

    for hdu in ext_mapping:
        for kw in keywords:
            try:
                del f[hdu].header[kw]
            except KeyError:
                pass
    ext_mapping = list(altwcs.mapFitsExt2HDUListInd(f, extname).values())
    ext_mapping.sort()
    for hdu in ext_mapping[::-1]:
        del f[hdu]
    f.close()
