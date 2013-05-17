from __future__ import division # confidence high
import os
import string

import numpy as np
import pywcs
import pyfits
from stsci.tools import fileutil as fu

altwcskw = ['WCSAXES', 'CRVAL', 'CRPIX', 'PC', 'CDELT', 'CD', 'CTYPE', 'CUNIT',
            'PV', 'PS']
altwcskw_extra = ['LATPOLE','LONPOLE','RESTWAV','RESTFRQ']

# file operations
def archiveWCS(fname, ext, wcskey=" ", wcsname=" ", reusekey=False):
    """
    Copy the primary WCS to the header as an alternate WCS
    with wcskey and name WCSNAME. It loops over all extensions in 'ext'

    Parameters
    ----------
    fname:  string or pyfits.HDUList
            a file name or a file object
    ext:    an int, a tuple, string, or list of integers or tuples (e.g.('sci',1))
            fits extensions to work with
            If a string is provided, it should specify the EXTNAME of extensions
            with WCSs to be archived
    wcskey: string "A"-"Z" or " "
            if " ": get next available key if wcsname is also " " or try
            to get a key from WCSNAME value
    wcsname: string
             Name of alternate WCS description
    reusekey: boolean
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

    if isinstance(fname, str):
        f = pyfits.open(fname, mode='update')
    else:
        f = fname

    if not _parpasscheck(f, ext, wcskey, wcsname):
        closefobj(fname, f)
        raise ValueError("Input parameters problem")

    # Interpret input 'ext' value to get list of extensions to process
    ext = _buildExtlist(f, ext)

    if not wcskey and not wcsname:
        raise KeyError("Either wcskey or wcsname should be specified")

    if wcsname.strip() == "":
        try:
            wcsname = readAltWCS(f, ext[0], wcskey=" ")['WCSNAME']
        except KeyError:
            pass
    wcsext = ext[0]
    if wcskey != " " and wcskey in wcskeys(f[wcsext].header) and not reusekey:
        closefobj(fname, f)
        raise KeyError("Wcskey %s is aready used. \
        Run archiveWCS() with reusekey=True to overwrite this alternate WCS. \
        Alternatively choose another wcskey with altwcs.available_wcskeys()." %wcskey)
    elif wcskey == " ":
        # wcsname exists, overwrite it if reuse is True or get the next key
        if wcsname.strip() in wcsnames(f[wcsext].header).values():
            if reusekey:
                # try getting the key from an existing WCS with WCSNAME
                wkey = getKeyFromName(f[wcsext].header, wcsname)
                wname = wcsname
                if wkey == ' ':
                    wkey = next_wcskey(f[wcsext].header)
                elif wkey is None:
                    closefobj(fname, f)
                    raise KeyError("Could not get a valid wcskey from wcsname %s" %wcsname)
            else:
                closefobj(fname, f)
                raise KeyError("Wcsname %s is aready used. \
                Run archiveWCS() with reusekey=True to overwrite this alternate WCS. \
                Alternatively choose another wcskey with altwcs.available_wcskeys() or\
                choose another wcsname." %wcsname)
        else:
            wkey = next_wcskey(f[wcsext].header)
            if wcsname.strip():
                wname = wcsname
            else:
                # determine which WCSNAME needs to be replicated in archived WCS
                wnames = wcsnames(f[wcsext].header)
                if 'O' in wnames: del wnames['O'] # we don't want OPUS/original
                if len(wnames) > 0:
                    if ' ' in wnames:
                        wname = wnames[' ']
                    else:
                        akeys = string.uppercase
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

    for e in ext:
        hwcs = readAltWCS(f,e,wcskey=' ')
        if hwcs is None:
            continue

        wcsnamekey = 'WCSNAME' + wkey
        f[e].header[wcsnamekey] = wname

        try:
            old_wcsname=hwcs.pop('WCSNAME')
        except:
            pass
        for k in hwcs.keys():

            key = k[:7] + wkey
            f[e].header[key] = hwcs[k]
    closefobj(fname, f)

def restore_from_to(f, fromext=None, toext=None, wcskey=" ", wcsname=" "):
    """
    Copy an alternate WCS from one extension as a primary WCS of another extension

    Reads in a WCS defined with wcskey and saves it as the primary WCS.
    Goes sequentially through the list of extensions in ext.
    Alternatively uses 'fromext' and 'toext'.


    Parameters
    ----------
    f:       string or pyfits.HDUList object
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
        fobj = pyfits.open(f, mode='update')
    else:
        fobj = f

    if not _parpasscheck(fobj, ext=None, wcskey=wcskey, fromext=fromext, toext=toext):
        closefobj(f, fobj)
        raise ValueError("Input parameters problem")

    # Interpret input 'ext' value to get list of extensions to process
    #ext = _buildExtlist(fobj, ext)

    if isinstance(toext, str):
        toext = [toext]

    # the case of an HDUList object in memory without an associated file

    #if fobj.filename() is not None:
    #        name = fobj.filename()

    simplefits = fu.isFits(fobj)[1] is 'simple'
    if simplefits:
        wcskeyext = 0
    else:
        wcskeyext = 1

    if wcskey == " ":
        if wcsname.strip():
            wkey = getKeyFromName(fobj[wcskeyext].header, wcsname)
            if not wkey:
                closefobj(f, fobj)
                raise KeyError("Could not get a key from wcsname %s ." % wcsname)
    else:
        if wcskey not in wcskeys(fobj, ext=wcskeyext):
            print "Could not find alternate WCS with key %s in this file" % wcskey
            closefobj(f, fobj)
            return
        wkey = wcskey

    countext = fu.countExtn(fobj, fromext)
    if not countext:
        raise KeyError("File does not have extension with extname %s", fromext)
    else:
        for i in range(1, countext+1):
            for toe in toext:
                _restore(fobj, fromextnum=i, fromextnam=fromext, toextnum=i, toextnam=toe, ukey=wkey)

    if fobj.filename() is not None:
        #fobj.writeto(name)
        closefobj(f, fobj)

def restoreWCS(f, ext, wcskey=" ", wcsname=" "):
    """
    Copy a WCS with key "WCSKEY" to the primary WCS

    Reads in a WCS defined with wcskey and saves it as the primary WCS.
    Goes sequentially through the list of extensions in ext.
    Alternatively uses 'fromext' and 'toext'.


    Parameters
    ----------
    f:       string or pyfits.HDUList object
             a file name or a file object
    ext:    an int, a tuple, string, or list of integers or tuples (e.g.('sci',1))
            fits extensions to work with
            If a string is provided, it should specify the EXTNAME of extensions
            with WCSs to be archived
    wcskey:  a charater
             "A"-"Z" - Used for one of 26 alternate WCS definitions.
             or " " - find a key from WCSNAMe value
    wcsname: string (optional)
             if given and wcskey is " ", will try to restore by WCSNAME value

    See Also
    --------
    archiveWCS - copy the primary WCS as an alternate WCS
    restore_from_to

    """
    if isinstance(f, str):
        fobj = pyfits.open(f, mode='update')
    else:
        fobj = f

    if not _parpasscheck(fobj, ext=ext, wcskey=wcskey):
        closefobj(f, fobj)
        raise ValueError("Input parameters problem")

    # Interpret input 'ext' value to get list of extensions to process
    ext = _buildExtlist(fobj, ext)


    # the case of an HDUList object in memory without an associated file

    #if fobj.filename() is not None:
    #        name = fobj.filename()

    simplefits = fu.isFits(fobj)[1] is 'simple'
    if simplefits:
        wcskeyext = 0
    else:
        wcskeyext = 1

    if wcskey == " ":
        if wcsname.strip():
            wkey = getKeyFromName(fobj[wcskeyext].header, wcsname)
            if not wkey:
                closefobj(f, fobj)
                raise KeyError("Could not get a key from wcsname %s ." % wcsname)
    else:
        if wcskey not in wcskeys(fobj, ext=wcskeyext):
            #print "Could not find alternate WCS with key %s in this file" % wcskey
            closefobj(f, fobj)
            return
        wkey = wcskey

    for e in ext:
        _restore(fobj, wkey, fromextnum=e, verbose=False)

    if fobj.filename() is not None:
        #fobj.writeto(name)
        closefobj(f, fobj)

def deleteWCS(fname, ext, wcskey=" ", wcsname=" "):
    """
    Delete an alternate WCS defined with wcskey.
    If wcskey is " " try to get a key from WCSNAME.

    Parameters
    ----------
    fname:   sting or a pyfits.HDUList object
    ext:    an int, a tuple, string, or list of integers or tuples (e.g.('sci',1))
            fits extensions to work with
            If a string is provided, it should specify the EXTNAME of extensions
            with WCSs to be archived
    wcskey:  one of 'A'-'Z' or " "
    wcsname: string
             Name of alternate WCS description
    """
    if isinstance(fname, str):
        fobj = pyfits.open(fname, mode='update')
    else:
        fobj = fname

    if not _parpasscheck(fobj, ext, wcskey, wcsname):
        closefobj(fname, fobj)
        raise ValueError("Input parameters problem")

    # Interpret input 'ext' value to get list of extensions to process
    ext = _buildExtlist(fobj, ext)
    # Do not allow deleting the original WCS.
    if wcskey == 'O':
        print "Wcskey 'O' is reserved for the original WCS and should not be deleted."
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
            raise KeyError("Could not get a key: wcsname '%s' not found in header." % wcsname)
    else:
        if wcskey not in wcskeys(fobj[wcskeyext].header):
            closefobj(fname, fobj)
            raise KeyError("Could not find alternate WCS with key %s in this file" % wcskey)
        wkey = wcskey

    prexts = []
    for i in ext:
        hdr = fobj[i].header
        hwcs = readAltWCS(fobj,i,wcskey=wkey)
        if hwcs is None:
            continue
        for k in hwcs.keys():
            del hdr[k]
            #del hdr['ORIENT'+wkey]
        prexts.append(i)
    if prexts != []:
        print 'Deleted all instances of WCS with key %s in extensions' % wkey, prexts
    else:
        print "Did not find WCS with key %s in any of the extensions" % wkey
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
    if not isinstance(ext,list):
        if isinstance(ext,str):
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
    if fromextnam:
        fromextension = (fromextnam, fromextnum)
    else:
        fromextension = fromextnum
    if toextnum:
        if toextnam:
            toextension = (toextnam, toextnum)
        else:
            toextension =toextnum
    else:
        toextension = fromextension

    hwcs = readAltWCS(fobj,fromextension,wcskey=ukey,verbose=verbose)
    if hwcs is None:
        return

    for k in hwcs.keys():
        key = k[:-1]
        if key in fobj[toextension].header.keys():
            fobj[toextension].header.update(key=key, value = hwcs[k])
        else:
            continue
    if key == 'O' and 'TDDALPHA' in fobj[toextension].header:
        fobj[toextension].header['TDDALPHA'] = 0.0
        fobj[toeztension].header['TDDBETA'] = 0.0
    if 'ORIENTAT' in fobj[toextension].header:
        norient = np.rad2deg(np.arctan2(hwcs['CD1_2'+'%s' %ukey], hwcs['CD2_2'+'%s' %ukey]))
        fobj[toextension].header.update(key='ORIENTAT', value=norient)

#header operations
def _check_headerpars(fobj, ext):
    if not isinstance(fobj, pyfits.Header) and not isinstance(fobj, pyfits.HDUList) \
                   and not isinstance(fobj, str):
        raise ValueError("Expected a file name, a file object or a header\n")

    if not isinstance(fobj, pyfits.Header):
        #raise ValueError("Expected a valid ext parameter when input is a file")
        if not isinstance(ext, int) and not isinstance(ext, tuple):
            raise ValueError("Expected ext to be a number or a tuple, e.g. ('SCI', 1)\n")

def _getheader(fobj, ext):
    if isinstance(fobj, str):
        hdr = pyfits.getheader(fobj,ext)
    elif isinstance(fobj, pyfits.Header):
        hdr = fobj
    else:
        hdr = fobj[ext].header
    return hdr

def readAltWCS(fobj, ext, wcskey=' ',verbose=False):
    """ Reads in alternate WCS from specified extension

    Parameters
    ----------
    fobj: string, pyfits.HDUList
        fits filename or pyfits file object
        containing alternate/primary WCS(s) to be converted
    wcskey: string [" ",A-Z]
        alternate/primary WCS key that will be replaced by the new key
    ext: int
        extension number

    Returns
    -------
    hdr: pyfits.Header
        header object with ONLY the keywords for specified alternate WCS
    """
    if isinstance(fobj, str):
        fobj = pyfits.open(fobj)

    hdr = _getheader(fobj,ext)
    try:
        nwcs = pywcs.WCS(hdr, fobj=fobj, key=wcskey)
    except KeyError:
        if verbose:
            print 'readAltWCS: Could not read WCS with key %s' %wcskey
            print '            Skipping %s[%s]' % (fobj.filename(), str(ext))
        return None
    hwcs = nwcs.to_header()

    if nwcs.wcs.has_cd():
        hwcs = pc2cd(hwcs, key=wcskey)

    return hwcs

def convertAltWCS(fobj,ext,oldkey=" ",newkey=' '):
    """
    Translates the alternate/primary WCS with one key to an alternate/primary WCS with
    another key.

    Parameters
    ----------
    fobj: string, pyfits.HDUList, or pyfits.Header
        fits filename, pyfits file object or pyfits header
        containing alternate/primary WCS(s) to be converted
    ext: int
        extension number
    oldkey: string [" ",A-Z]
        alternate/primary WCS key that will be replaced by the new key
    newkey: string [" ",A-Z]
        new alternate/primary WCS key

    Returns
    -------
    hdr: pyfits.Header
        header object with keywords renamed from oldkey to newkey
    """
    hdr = readAltWCS(fobj,ext,wcskey=oldkey)
    if hdr is None:
        return None
    # Converting WCS to new key
    for card in hdr:
        if oldkey == ' ' or oldkey == '':
            cname = card
        else:
            cname = card.rstrip(oldkey)
        hdr.rename_key(card,cname+newkey,force=True)

    return hdr

def wcskeys(fobj, ext=None):
    """
    Returns a list of characters used in the header for alternate
    WCS description with WCSNAME keyword

    Parameters
    ----------
    fobj: string, pyfits.HDUList or pyfits.Header
          fits file name, pyfits file object or pyfits header
    ext: int or None
         extension number
         if None, fobj must be a header
    """
    _check_headerpars(fobj, ext)
    hdr = _getheader(fobj, ext)
    names = hdr["WCSNAME*"]
    d = []
    for key in names:
        wkey = key.replace('WCSNAME','')
        if wkey == '': wkey = ' '
        d.append(wkey)
    return d

def wcsnames(fobj, ext=None):
    """
    Returns a dictionary of wcskey: WCSNAME pairs

    Parameters
    ----------
    fobj: string, pyfits.HDUList or pyfits.Header
          fits file name, pyfits file object or pyfits header
    ext: int or None
         extension number
         if None, fobj must be a header

    """
    _check_headerpars(fobj, ext)
    hdr = _getheader(fobj, ext)
    names = hdr["WCSNAME*"]
    d = {}
    for keyword, value in names.items():
        wkey = keyword.replace('WCSNAME','')
        if wkey == '': wkey = ' '
        d[wkey] = value
    return d

def available_wcskeys(fobj, ext=None):
    """
    Returns a list of characters which are not used in the header
    with WCSNAME keyword. Any of them can be used to save a new
    WCS.

    Parameters
    ----------
    fobj: string, pyfits.HDUList or pyfits.Header
          fits file name, pyfits file object or pyfits header
    ext: int or None
         extension number
         if None, fobj must be a header
    """
    _check_headerpars(fobj, ext)
    hdr = _getheader(fobj, ext)
    all_keys = list(string.ascii_uppercase)
    used_keys = wcskeys(hdr)
    try:
        used_keys.remove(" ")
    except ValueError:
        pass
    [all_keys.remove(key) for key in used_keys]
    return all_keys

def next_wcskey(fobj, ext=None):
    """
    Returns next available character to be used for an alternate WCS

    Parameters
    ----------
    fobj: string, pyfits.HDUList or pyfits.Header
          fits file name, pyfits file object or pyfits header
    ext: int or None
         extension number
         if None, fobj must be a header
    """
    _check_headerpars(fobj, ext)
    hdr = _getheader(fobj, ext)
    allkeys = available_wcskeys(hdr)
    if allkeys != []:
        return allkeys[0]
    else:
        return None

def getKeyFromName(header, wcsname):
    """
    If WCSNAME is found in header, return its key, else return
    None. This is used to update an alternate WCS
    repeatedly and not generate new keys every time.

    Parameters
    ----------
    header:  pyfits.Header
    wcsname: str
             Value of WCSNAME
    """
    wkey = None
    names = wcsnames(header)
    for item in names.items():
        if item[1].lower() == wcsname.lower():
            wkey = item[0]
            break
    return wkey

def pc2cd(hdr, key=' '):
    """
    Convert a CD PC matrix to a CD matrix.

    WCSLIB (and PyWCS) recognizes CD keywords as input
    but converts them and works internally with the PC matrix.
    to_header() returns the PC matrix even if the i nput was a
    CD matrix. To keep input and output consistent we check
    for has_cd and convert the PC back to CD.

    Parameters
    ----------
    hdr: pyfits.Header

    """
    for c in ['1_1', '1_2', '2_1', '2_2']:
        try:
            val = hdr['PC'+c+'%s' % key]
            del hdr['PC'+c+ '%s' % key]
        except KeyError:
            if c=='1_1' or c == '2_2':
                val = 1.
            else:
                val = 0.
        hdr.update(key='CD'+c+'%s' %key, value=val)
    return hdr

def _parpasscheck(fobj, ext, wcskey, fromext=None, toext=None, reusekey=False):
    """
    Check input parameters to altwcs functions

    fobj:    string or pyfits.HDUList object
             a file name or a file object
    ext:     an int, a tuple, a python list of integers or a python list
             of tuples (e.g.('sci',1))
             fits extensions to work with
    wcskey:  a charater
             "A"-"Z" or " "- Used for one of 26 alternate WCS definitions
    wcsname: string (optional)
             if given and wcskey is " ", will try to restore by WCSNAME value
    reusekey: boolean
             A flag which indicates whether to reuse a wcskey in the header
    """
    if not isinstance(fobj,pyfits.HDUList):
        print "First parameter must be a fits file object or a file name."
        return False

    # first one covers the case of an object created in memory
    # (e.g. headerlet) for which fileinfo returns None
    if fobj.fileinfo(0) is None:
        pass
    else:
        # an HDUList object with associated file
        if fobj.fileinfo(0)['filemode'] is not 'update':
            print "First parameter must be a file name or a file object opened in 'update' mode."
            return False

    if not isinstance(ext, int) and not isinstance(ext, tuple) \
        and not isinstance(ext,str) \
        and not isinstance(ext, list) and ext is not None:
        print "Ext must be integer, tuple, string,a list of int extension numbers, \n\
        or a list of tuples representing a fits extension, for example ('sci', 1)."
        return False

    if not isinstance(fromext, str) and fromext is not None:
        print "fromext must be a string representing a valid extname"
        return False

    if not isinstance(toext, list) and not isinstance(toext, str) and \
                     toext is not None :
        print "toext must be a string or a list of strings representing extname"
        return False

    if len(wcskey) != 1:
        print 'Parameter wcskey must be a character - one of "A"-"Z" or " "'
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

    if not isinstance(fname, pyfits.HDUList):
        f = pyfits.open(fname)
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
