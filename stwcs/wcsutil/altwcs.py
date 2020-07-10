import string

import numpy as np
from astropy import wcs as pywcs
from astropy.io import fits
from stsci.tools import fileutil as fu

from astropy import log
from astropy.utils.decorators import deprecated


default_log_level = log.getEffectiveLevel()


__all__ = ["archiveWCS", "available_wcskeys", "convertAltWCS", "deleteWCS", "next_wcskey",
           "pc2cd", "readAltWCS", "restoreWCS", "wcskeys", "wcsnames",
           "wcs_from_key"]


altwcskw = ['WCSAXES', 'CRVAL', 'CRPIX', 'PC', 'CDELT', 'CD', 'CTYPE', 'CUNIT',
            'PV', 'PS']
altwcskw_extra = ['LATPOLE', 'LONPOLE', 'RESTWAV', 'RESTFRQ']

# List non-standard WCS keywords (such as those created, e.g., by TweakReg)
# that need to be archived/restored with the rest of WCS here:
STWCS_KWDS = ['WCSTYPE', 'RMS_RA', 'RMS_DEC', 'NMATCH', 'FITNAME']

# file operations


def archiveWCS(fname, ext, wcskey=" ", wcsname=" ", reusekey=False):
    """
    Copy the primary WCS to the header as an alternate WCS
    with wcskey and name WCSNAME. It loops over all extensions in 'ext'

    Parameters
    ----------
    fname :  string or `astropy.io.fits.HDUList`
        file name or a file object
    ext :    int, tuple, str, or list of integers or tuples (e.g.('sci',1))
        fits extensions to work with
        If a string is provided, it should specify the EXTNAME of extensions
        with WCSs to be archived
    wcskey : string "A"-"Z" or " "
        if ``' '``: get next available key if wcsname is also ``' '`` or try
        to get a key from WCSNAME value
    wcsname : string
        Name of alternate WCS description
    reusekey : boolean
        if True - overwrites a WCS with the same key

    Examples
    --------
    Copy the primary WCS of an in memory headrlet object to an
    alternate WCS with key 'T'

    >>> hlet=headerlet.createHeaderlet('junk.fits', 'hdr1.fits')
    >>> altwcs.wcskeys(hlet[1].header)
    ['A']
    >>> altwcs.archiveWCS(hlet, ext=[('SIPWCS',1),('SIPWCS',2)], wcskey='T')
    >>> altwcs.wcskeys(hlet[1].header)
    ['A', 'T']


    See Also
    --------
    wcsutil.restoreWCS: Copy an alternate WCS to the primary WCS

    """
    if wcsname.strip() == '':
        wcsname = ' '

    if len(wcskey) != 1 or wcskey.strip() not in string.ascii_letters:
        raise ValueError(
            "Parameter 'wcskey' must be a character - one of 'A'-'Z' or ' '."
        )

    if isinstance(fname, str):
        f = fits.open(fname, mode='update')
    else:
        f = fname

    if not _parpasscheck(f, ext, wcskey, wcsname):
        closefobj(fname, f)
        raise ValueError("Input parameters problem")

    # Interpret input 'ext' value to get list of extensions to process
    ext = _buildExtlist(f, ext)

    if wcsname == ' ':
        try:
            wcsname = wcs_from_key(f, ext[0])['WCSNAME']
        except KeyError:
            pass

    wcsext = ext[0]
    if wcskey != " " and wcskey in wcskeys(f[wcsext].header) and not reusekey:
        closefobj(fname, f)
        raise KeyError(f"Wcskey '{wcskey}' is aready used. \
        Run archiveWCS() with reusekey=True to overwrite this alternate WCS. \
        Alternatively choose another wcskey with altwcs.available_wcskeys().")

    elif wcskey == " ":
        # wcsname exists, overwrite it if reuse is True or get the next key
        if wcsname in _alt_wcs_names(f[wcsext].header).values():
            if reusekey:
                # try getting the key from an existing WCS with WCSNAME
                wkey = getKeyFromName(f[wcsext].header, wcsname)
                wname = wcsname
                if wkey == ' ':
                    wkey = _next_wcskey(f[wcsext].header)
                elif wkey is None:
                    closefobj(fname, f)
                    raise KeyError(f"Could not get a valid wcskey from wcsname '{wcsname:s}'")
            else:
                closefobj(fname, f)
                raise KeyError(
                    f"Wcsname '{wcsname}' is aready used. Run archiveWCS() "
                    "with reusekey=True to overwrite this alternate WCS. "
                    "Alternatively choose another wcskey with "
                    "altwcs.available_wcskeys() or choose another wcsname."
                )

        else:
            wkey = _next_wcskey(f[wcsext].header)
            if wcsname != ' ':
                wname = wcsname
            else:
                # determine which WCSNAME needs to be replicated in archived WCS
                wnames = _alt_wcs_names(f[wcsext].header)

                if len(wnames) > 0:
                    if ' ' in wnames:
                        wname = wnames[' ']
                    else:
                        akeys = string.ascii_uppercase
                        wname = "DEFAULT"
                        for key in akeys[-1::]:
                            if key in wnames:
                                wname = wnames
                                break
                else:
                    wname = "DEFAULT"
    else:
        wkey = wcskey
        wname = wcsname

    log.setLevel('WARNING')

    for e in ext:
        hwcs = wcs_from_key(f, e, from_key=' ', to_key=wkey)

        if hwcs is None:
            continue

        hwcs.set(f'WCSNAME{wkey:.1s}', wname, 'Coordinate system title', before=0)

        f[e].header.update(hwcs)

    log.setLevel(default_log_level)
    closefobj(fname, f)


def restore_from_to(f, fromext=None, toext=None, wcskey=" ", wcsname=" "):
    """
    Copy an alternate WCS from one extension as a primary WCS of another extension

    Reads in a WCS defined with wcskey and saves it as the primary WCS.
    Goes sequentially through the list of extensions in ext.
    Alternatively uses 'fromext' and 'toext'.


    Parameters
    ----------
    f:       string or `astropy.io.fits.HDUList`
             a file name or a file object
    fromext: string
             extname from which to read in the alternate WCS, for example 'SCI'
    toext:   string or python list
             extname or a list of extnames to which the WCS will be copied as
             primary, for example ['SCI', 'ERR', 'DQ']
    wcskey:  a charater
             "A"-"Z" - Used for one of 26 alternate WCS definitions.
             or " " - find a key from WCSNAMe value
    wcsname: string (optional)
             if given and wcskey is " ", will try to restore by WCSNAME value

    See Also
    --------
    archiveWCS - copy the primary WCS as an alternate WCS
    restoreWCS - Copy a WCS with key "WCSKEY" to the primary WCS

    """
    if isinstance(f, str):
        fobj = fits.open(f, mode='update')
    else:
        fobj = f

    if not _parpasscheck(fobj, ext=None, wcskey=wcskey, fromext=fromext, toext=toext):
        closefobj(f, fobj)
        raise ValueError("Input parameters problem")

    # Interpret input 'ext' value to get list of extensions to process
    # ext = _buildExtlist(fobj, ext)

    if isinstance(toext, str):
        toext = [toext]

    # the case of an HDUList object in memory without an associated file

    # if fobj.filename() is not None:
    #        name = fobj.filename()

    wcskeyext = 0 if fu.isFits(fobj)[1] == 'simple' else 1

    if wcskey == " ":
        if wcsname.strip():
            wkey = getKeyFromName(fobj[wcskeyext].header, wcsname)
            if not wkey:
                closefobj(f, fobj)
                raise KeyError(f"Could not get a key from wcsname '{wcsname}'.")
    else:
        if wcskey not in wcskeys(fobj, ext=wcskeyext):
            print(f"Could not find alternate WCS with key '{wcskey}' in this file")
            closefobj(f, fobj)
            return
        wkey = wcskey

    countext = fu.countExtn(fobj, fromext)
    if countext:
        for i in range(1, countext + 1):
            for toe in toext:
                _restore(fobj, fromextnum=i, fromextnam=fromext, toextnum=i,
                         toextnam=toe, ukey=wkey)
    else:
        raise KeyError(f"File does not have extension with extname {fromext:s}")

    if fobj.filename() is not None:
        closefobj(f, fobj)


def restoreWCS(f, ext, wcskey=" ", wcsname=" "):
    """
    Copy a WCS with key "WCSKEY" to the primary WCS

    Reads in a WCS defined with wcskey and saves it as the primary WCS.
    Goes sequentially through the list of extensions in ext.
    Alternatively uses 'fromext' and 'toext'.


    Parameters
    ----------
    f : str or `astropy.io.fits.HDUList`
        file name or a file object
    ext : int, tuple, str, or list of integers or tuples (e.g.('sci',1))
        fits extensions to work with
        If a string is provided, it should specify the EXTNAME of extensions
        with WCSs to be archived
    wcskey : str
        "A"-"Z" - Used for one of 26 alternate WCS definitions.
        or " " - find a key from WCSNAMe value
    wcsname : str
        (optional) if given and wcskey is " ", will try to restore by WCSNAME value

    See Also
    --------
    archiveWCS - copy the primary WCS as an alternate WCS
    restore_from_to

    """
    if isinstance(f, str):
        fobj = fits.open(f, mode='update')
    else:
        fobj = f

    if not _parpasscheck(fobj, ext=ext, wcskey=wcskey):
        closefobj(f, fobj)
        raise ValueError("Input parameters problem")

    # Interpret input 'ext' value to get list of extensions to process
    ext = _buildExtlist(fobj, ext)

    # the case of an HDUList object in memory without an associated file

    wcskeyext = 0 if fu.isFits(fobj)[1] == 'simple' else 1

    if wcskey == " ":
        if wcsname.strip():
            wcskey = getKeyFromName(fobj[wcskeyext].header, wcsname)
            if not wcskey:
                closefobj(f, fobj)
                raise KeyError(f"Could not get a key from wcsname '{wcsname}'.")

    for e in ext:
        if wcskey in wcskeys(fobj, ext=e):
            _restore(fobj, wcskey, fromextnum=e, verbose=False)

    if fobj.filename() is not None:
        closefobj(f, fobj)


def deleteWCS(fname, ext, wcskey=" ", wcsname=" "):
    """
    Delete an alternate WCS defined with wcskey.
    If wcskey is " " try to get a key from WCSNAME.

    Parameters
    ----------
    fname : str or a `astropy.io.fits.HDUList`
    ext : int, tuple, str, or list of integers or tuples (e.g.('sci',1))
        fits extensions to work with
        If a string is provided, it should specify the EXTNAME of extensions
        with WCSs to be archived
    wcskey : str
        one of 'A'-'Z' or " "
    wcsname : str
        Name of alternate WCS description
    """
    if isinstance(fname, str):
        fobj = fits.open(fname, mode='update')
    else:
        fobj = fname

    if not _parpasscheck(fobj, ext, wcskey, wcsname):
        closefobj(fname, fobj)
        raise ValueError("Input parameters problem")

    # Interpret input 'ext' value to get list of extensions to process
    ext = _buildExtlist(fobj, ext)
    # Do not allow deleting the original WCS.
    if wcskey == 'O':
        print("Wcskey 'O' is reserved for the original WCS and should not be deleted.")
        closefobj(fname, fobj)
        return

    wcskeyext = ext[0]

    if not wcskeys and not wcsname:
        raise KeyError("Either wcskey or wcsname should be specified")

    if wcskey == " ":
        # try getting the key from WCSNAME
        wkey = getKeyFromName(fobj[wcskeyext].header, wcsname)
        if not wkey:
            closefobj(fname, fobj)
            raise KeyError(f"Could not get a key: wcsname '{wcsname}' not found in header.")
    else:
        if wcskey not in wcskeys(fobj[wcskeyext].header):
            closefobj(fname, fobj)
            raise KeyError(f"Could not find alternate WCS with key '{wcskey}' in this file")
        wkey = wcskey

    prexts = []
    for i in ext:
        hdr = fobj[i].header
        hwcs = wcs_from_key(fobj, i, from_key=wkey)
        if hwcs is not None:
            for k in hwcs:
                if k in hdr:
                    del hdr[k]
            prexts.append(i)

    if prexts:
        print(f'Deleted all instances of WCS with key {wkey:s} in extensions {prexts}')
    else:
        print(f"Did not find WCS with key {wkey:s} in any of the extensions {prexts}")
    closefobj(fname, fobj)


def _buildExtlist(fobj, ext):
    """
    Utility function to interpret the provided value of 'ext' and return a list
    of 'valid' values which can then be used by the rest of the functions in
    this module.

    Parameters
    ----------
    fobj: HDUList
        file to be examined
    ext:    an int, a tuple, string, list of integers or tuples (e.g.('sci',1))
            fits extensions to work with
            If a string is provided, it should specify the EXTNAME of extensions
            with WCSs to be archived
    """
    if not isinstance(ext, list):
        if isinstance(ext, str):
            extstr = ext
            ext = []
            for extn in range(1, len(fobj)):
                if 'extname' in fobj[extn].header and fobj[extn].header['extname'] == extstr:
                    ext.append(extn)
        elif isinstance(ext, int) or isinstance(ext, tuple):
            ext = [ext]
        else:
            raise KeyError("Valid extensions in 'ext' parameter need to be specified.")
    return ext


def _restore(fobj, ukey, fromextnum,
             toextnum=None, fromextnam=None, toextnam=None, verbose=True):
    """
    fobj: string of HDUList
    ukey: string 'A'-'Z'
          wcs key
    fromextnum: int
                extver of extension from which to copy WCS
    fromextnam: string
                extname of extension from which to copy WCS
    toextnum: int
              extver of extension to which to copy WCS
    toextnam: string
              extname of extension to which to copy WCS
    """
    # create an extension tuple, e.g. ('SCI', 2)
    fromextension = (fromextnam, fromextnum) if fromextnam else fromextnum

    if toextnum:
        if toextnam:
            toextension = (toextnam, toextnum)
        else:
            toextension = toextnum
    else:
        toextension = fromextension

    hwcs = wcs_from_key(fobj, fromextension, from_key=ukey, to_key=' ')
    fobj[toextension].header.update(hwcs)

    if ukey == 'O' and 'TDDALPHA' in fobj[toextension].header:
        fobj[toextension].header['TDDALPHA'] = 0.0
        fobj[toextension].header['TDDBETA'] = 0.0

    if 'ORIENTAT' in fobj[toextension].header:
        norient = np.rad2deg(np.arctan2(hwcs[f'CD1_2'], hwcs[f'CD2_2']))
        fobj[toextension].header['ORIENTAT'] = norient

    # Reset 2014 TDD keywords prior to computing new values (if any are computed)
    for kw in ['TDD_CYA', 'TDD_CYB', 'TDD_CXA', 'TDD_CXB']:
        if kw in fobj[toextension].header:
            fobj[toextension].header[kw] = 0.0


# header operations


def _check_headerpars(fobj, ext):
    if not isinstance(fobj, fits.Header) and not isinstance(fobj, fits.HDUList) \
            and not isinstance(fobj, str):
        raise ValueError("Expected a file name, a file object or a header\n")

    if not isinstance(fobj, fits.Header):
        if not isinstance(ext, int) and not isinstance(ext, tuple):
            raise ValueError("Expected ext to be a number or a tuple, e.g. ('SCI', 1)\n")


def _getheader(fobj, ext):
    if isinstance(fobj, str):
        hdr = fits.getheader(fobj, ext)
    elif isinstance(fobj, fits.Header):
        hdr = fobj
    else:
        hdr = fobj[ext].header
    return hdr


def wcs_from_key(fobj, ext, from_key=' ', to_key=None):
    """
    Read in WCS with a given ``from_key`` from the specified extension of
    the input file object and return a FITS header representing this WCS
    using desired WCS key specified by ``to_key``.

    Parameters
    ----------
    fobj : str, `astropy.io.fits.HDUList`
        FITS filename or `astropy.io.fits.HDUList` object containing
        a header with an alternate/primary WCS to be read.

    ext : int, str or tuple of (str, int)
        Extension specification identifying image HDU from which WCS should be
        loaded. If ``ext`` is a tuple, it is of the form ``(extname, extver)``
        where ``extname`` is a `str` extension name and ``extver`` is
        an integer extension version.
        An integer ``ext`` indicates "extension number". Finally, a single
        `str` extension name is interpretted as ``(ext, 1)``.

    from_key : {' ', 'A'-'Z'}
        A 1 character string that is either a space character indicating the
        primary WCS, or one of the 26 ASCII letters (``'A'``-``'Z'``)
        indicating alternate WCS to be loaded from specified header.

    to_key : {None, ' ', 'A'-'Z'}
        The key of the primary/alternate WCS to be used in the returned header.
        When ``to_key`` is `None`, the returned header describes a WCS with the
        same key as the one read in using ``from_key``. A space character or
        a single ASCII letter indicates the key to be used for the returned
        WCS (see ``from_key`` for details).

    Returns
    -------
    hdr: astropy.io.fits.Header
        Header object with FITS representation for specified primary or
        alternate WCS.

    """
    if len(from_key) != 1 or from_key.strip() not in string.ascii_uppercase:
        raise ValueError(
            "Parameter 'from_key' must be a character - one of 'A'-'Z' or ' '."
        )

    if to_key is None:
        to_key = from_key

    elif len(to_key) != 1 or to_key.strip() not in string.ascii_uppercase:
        raise ValueError(
            "Parameter 'to_key' must be a character - one of 'A'-'Z' or ' '."
        )

    if isinstance(fobj, str):
        fobj = fits.open(fobj)
        close_fobj = True
    else:
        close_fobj = False

    hdr = _getheader(fobj, ext)

    try:
        w = pywcs.WCS(hdr, fobj=fobj, key=from_key)
    except KeyError:
        log.warning(f'wcs_from_key: Could not read WCS with key {from_key:s}')
        log.warning(f'              Skipping {fobj.filename():s}[{ext}]')
        return fits.Header()
    finally:
        if close_fobj:
            fobj.close()

    hwcs = w.to_header(key=to_key)

    if w.wcs.has_cd():
        hwcs = pc2cd(hwcs, key=to_key)

    # preserve CTYPE values:
    for k, ctype in enumerate(w.wcs.ctype):
        hwcs[f'CTYPE{k + 1:d}{to_key:.1s}'] = ctype

    # include non-standard (i.e., tweakreg-specific) keywords
    from_key_s = from_key.strip()
    to_key_s = to_key.strip()
    for kwd in STWCS_KWDS:
        from_kwd = kwd + from_key_s
        if from_kwd in hdr:
            # preserves comments and empty/null str values
            idx = hdr.index(from_kwd)
            hdri = hdr[idx:idx+1]
            hdri.rename_keyword(from_kwd, kwd + to_key_s, force=True)
            hwcs.update(hdri)

    return hwcs


@deprecated(since='1.5.4', message='', name='readAltWCS',
            alternative='wcs_from_key')
def readAltWCS(fobj, ext, wcskey=' ', verbose=False):
    """
    Reads in alternate primary WCS from specified extension.

    Parameters
    ----------
    fobj : str, `astropy.io.fits.HDUList`
        fits filename or fits file object
        containing alternate/primary WCS(s) to be converted
    wcskey : str
        [" ",A-Z]
        alternate/primary WCS key that will be replaced by the new key
    ext : int
        fits extension number

    Returns
    -------
    hdr: fits.Header
        header object with ONLY the keywords for specified alternate WCS
    """
    return wcs_from_key(fobj, ext, from_key=wcskey)


@deprecated(since='1.5.4', message='', name='convertAltWCS',
            alternative='wcs_from_key')
def convertAltWCS(fobj, ext, oldkey=' ', newkey=' '):
    """
    Translates the alternate/primary WCS with one key to an alternate/primary WCS with
    another key.

    Parameters
    ----------
    fobj : str, `astropy.io.fits.HDUList`, or `astropy.io.fits.Header`
        fits filename, fits file object or fits header
        containing alternate/primary WCS(s) to be converted
    ext : int
        extension number
    oldkey : str
        [" ",A-Z]
        alternate/primary WCS key that will be replaced by the new key
    newkey : str
        [" ",A-Z]
        new alternate/primary WCS key

    Returns
    -------
    hdr: `astropy.io.fits.Header`
        header object with keywords renamed from oldkey to newkey
    """
    return wcs_from_key(fobj, ext, from_key=wcskey, to_key=newkey)


def wcskeys(fobj, ext=None):
    """
    Returns a list of characters used in the header for alternate
    WCS description with WCSNAME keyword

    Parameters
    ----------
    fobj : str, `astropy.io.fits.HDUList` or `astropy.io.fits.Header`
          fits file name, fits file object or fits header
    ext : int or None
         extension number
         if None, fobj must be a header
    """
    _check_headerpars(fobj, ext)
    hdr = _getheader(fobj, ext)

    wcs_kwd_list = ['WCSNAME', 'CTYPE1', 'CRPIX1', 'CRVAL1', 'RADESYS', 'LONPOLE']

    wkeys = set()

    # check primary:
    for kwd in wcs_kwd_list:
        if kwd in hdr:
            wkeys.add(' ')
            break

    # check Alt WCS:
    for kwd in wcs_kwd_list:
        alt_kwds = hdr[kwd + '?']
        alt_keys = [key[-1].upper() for key in alt_kwds if key[-1] in string.ascii_letters]
        wkeys.update(alt_keys)

    return sorted(wkeys)


def _alt_wcs_names(hdr, del_opus=True):
    """ Return a dictionary of all alternate WCS keys with names except ``OPUS``

    Parameters
    ----------
    hdr : astropy.io.fits.Header
        An image header.

    del_opus : bool
        Indicates whether to remove ``OPUS`` entry (WCS key ``'O'``) from
        returned key-name dictionary.

    Returns
    -------
    wnames : dict
        A dictionary of **Alt** WCS keys and names (as values):

    """
    names = hdr["WCSNAME?"]

    if del_opus and 'WCSNAMEO' in names and names['WCSNAMEO'] == 'OPUS':
        # remove OPUS name
        del names['WCSNAMEO']

    wnames = {kwd[-1].upper(): val for kwd, val in names.items()
              if kwd[-1] in string.ascii_letters}

    return wnames


def wcsnames(fobj, ext=None):
    """
    Returns a dictionary of wcskey: WCSNAME pairs

    Parameters
    ----------
    fobj : stri, `astropy.io.fits.HDUList` or `astropy.io.fits.Header`
          fits file name, fits file object or fits header
    ext : int or None
         extension number
         if None, fobj must be a header

    """
    _check_headerpars(fobj, ext)
    hdr = _getheader(fobj, ext)
    return _alt_wcs_names(hdr, del_opus=False)


def available_wcskeys(fobj, ext=None):
    """
    Returns a list of characters which are not used in the header
    with WCSNAME keyword. Any of them can be used to save a new
    WCS.

    Parameters
    ----------
    fobj : str, `astropy.io.fits.HDUList` or `astropy.io.fits.Header`
        fits file name, fits file object or fits header
    ext : int or None
        extension number
        if None, fobj must be a header
    """
    _check_headerpars(fobj, ext)
    hdr = _getheader(fobj, ext)
    return sorted(set(string.ascii_uppercase).difference(wcskeys(hdr)))


def next_wcskey(fobj, ext=None):
    """
    Returns next available character to be used for an alternate WCS

    Parameters
    ----------
    fobj : str, `astropy.io.fits.HDUList` or `astropy.io.fits.Header`
        fits file name, fits file object or fits header
    ext : int or None
        extension number
        if None, fobj must be a header
    """
    _check_headerpars(fobj, ext)
    hdr = _getheader(fobj, ext)
    return _next_wcskey(hdr)


def _next_wcskey(hdr):
    """
    Returns next available character to be used for an alternate WCS

    Parameters
    ----------
    hdr : `astropy.io.fits.Header`
        FITS header.

    Returns
    -------
    key : str, None
        Next available character to be used as an alternate WCS key or `None`
        if none is available.

    """
    all_keys = available_wcskeys(hdr)
    key = all_keys[0] if all_keys else None
    return key


def getKeyFromName(header, wcsname):
    """
    If WCSNAME is found in header, return its key, else return
    None. This is used to update an alternate WCS
    repeatedly and not generate new keys every time.

    Parameters
    ----------
    header :  `astropy.io.fits.Header`
    wcsname : str
        value of WCSNAME
    """
    names = wcsnames(header)
    wkeys = [key for key, name in names.items() if name.lower() == wcsname.lower()]
    if wkeys:
        wkey = max(wkeys)
    else:
        wkey = None
    return wkey


def pc2cd(hdr, key=' '):
    """
    Convert a PC matrix to a CD matrix.

    WCSLIB (and PyWCS) recognizes CD keywords as input
    but converts them and works internally with the PC matrix.
    to_header() returns the PC matrix even if the input was a
    CD matrix. To keep input and output consistent we check
    for has_cd and convert the PC back to CD.

    Parameters
    ----------
    hdr: `astropy.io.fits.Header`

    """
    key = key.strip()
    cdelt1 = hdr.pop(f'CDELT1{key:.1s}', 1)
    cdelt2 = hdr.pop(f'CDELT2{key:.1s}', 1)
    hdr[f'CD1_1{key:.1s}'] = (cdelt1 * hdr.pop(f'PC1_1{key:.1s}', 1),
                              'partial of first axis coordinate w.r.t. x')
    hdr[f'CD1_2{key:.1s}'] = (cdelt1 * hdr.pop(f'PC1_2{key:.1s}', 0),
                              'partial of first axis coordinate w.r.t. y')
    hdr[f'CD2_1{key:.1s}'] = (cdelt2 * hdr.pop(f'PC2_1{key:.1s}', 0),
                              'partial of second axis coordinate w.r.t. x')
    hdr[f'CD2_2{key:.1s}'] = (cdelt2 * hdr.pop(f'PC2_2{key:.1s}', 1),
                              'partial of second axis coordinate w.r.t. y')
    return hdr


def _parpasscheck(fobj, ext, wcskey, fromext=None, toext=None, reusekey=False):
    """
    Check input parameters to altwcs functions

    fobj : str or `astropy.io.fits.HDUList` object
        a file name or a file object
    ext : int, a tuple,a python list of integers or a python list
        of tuples (e.g.('sci',1))
        fits extensions to work with
    wcskey : str
        "A"-"Z" or " "- Used for one of 26 alternate WCS definitions
    wcsname : str
        (optional)
        if given and wcskey is " ", will try to restore by WCSNAME value
    reusekey : bool
        A flag which indicates whether to reuse a wcskey in the header
    """
    if not isinstance(fobj, fits.HDUList):
        print("First parameter must be a fits file object or a file name.")
        return False

    # first one covers the case of an object created in memory
    # (e.g. headerlet) for which fileinfo returns None
    if fobj.fileinfo(0) is None:
        pass
    else:
        # an HDUList object with associated file
        if fobj.fileinfo(0)['filemode'] != 'update':
            print("First parameter must be a file name or a file object opened in 'update' mode.")
            return False

    if not isinstance(ext, int) and not isinstance(ext, tuple) \
        and not isinstance(ext, str) \
        and not isinstance(ext, list) and ext is not None:
        print("Ext must be integer, tuple, string,a list of int extension "
              "numbers,\nor a list of tuples representing a fits extension, "
              "for example ('sci', 1).")
        return False

    if not isinstance(fromext, str) and fromext is not None:
        print("fromext must be a string representing a valid extname")
        return False

    if not isinstance(toext, list) and not isinstance(toext, str) and \
                        toext is not None:
        print("toext must be a string or a list of strings representing extname")
        return False

    if len(wcskey) > 1 or wcskey.strip() not in string.ascii_letters:
        print('Parameter wcskey must be a character - one of "A"-"Z" or " "')
        return False

    return True


def closefobj(fname, f):
    """
    Functions in this module accept as input a file name or a file object.
    If the input was a file name (string) we close the object. If the user
    passed a file object we leave it to the user to close it.
    """
    if isinstance(fname, str):
        f.close()


def mapFitsExt2HDUListInd(fname, extname):
    """
    Map FITS extensions with 'EXTNAME' to HDUList indexes.
    """

    if not isinstance(fname, fits.HDUList):
        f = fits.open(fname)
        close_file = True
    else:
        f = fname
        close_file = False

    d = {}
    for hdu in f:
        if 'EXTNAME' in hdu.header and hdu.header['EXTNAME'] == extname:
            extver = hdu.header['EXTVER']
            d[(extname, extver)] = f.index_of((extname, extver))

    if close_file:
        f.close()
    return d
