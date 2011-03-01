from __future__ import division # confidence high

import os.path, string
import pywcs
import numpy as np
import pyfits

altwcskw = ['WCSAXES', 'CRVAL', 'CRPIX', 'PC', 'CDELT', 'CD', 'CTYPE', 'CUNIT', 'PV', 'PS']
# file operations
def archiveWCS(fname, ext, wcskey=" ", wcsname=" ", clobber=False):
    """
    Copy the primary WCS to the hader as an alternate WCS 
    with wcskey and name WCSNAME. It loops over all extensions in 'ext'   
    
    Parameters
    ----------
    fname:  string or pyfits.HDUList
            a file name or a file object
    ext:    an int, a tuple, a python list of integers or a python list of tuples (e.g.('sci',1))
            fits extensions to work with
    wcskey: string "A"-"Z" or " " 
            if " ": get next available key if wcsname is also " " or try 
            to get a key from WCSNAME value 
    wcsname: string
             Name of alternate WCS description
    clobber: boolean
             if Ture - overwrites a WCS with the same key
             
    See Also
    --------
    wcsutils.restoreWCS: Copy an alternate WCS to the primary WCS
    
    """
    
    if isinstance(fname, str):
        f = pyfits.open(fname, mode='update')
    else:
        f = fname

    if not parpasscheck(f, ext, wcskey, wcsname):
        closefobj(fname, f)
        return
    
    if isinstance(ext, int) or isinstance(ext, tuple):
        ext = [ext]
        
    if wcskey == " ": 
        # try getting the key from WCSNAME
        if not wcsname.strip():
            wkey = next_wcskey(f[1].header)
        else:
            wkey = getKeyFromName(f[1].header, wcsname)
    else:
        if wcskey not in available_wcskeys(f[1].header):
            # reuse wcsname
            if not wcsname.strip(): 
                wcsname = f[1].header["WCSNAME"+wcskey]
                wname = wcsname
            else:
                wkey = wcskey
                wname = wcsname
        else:
            wkey = wcskey 
            wname = wcsname
            
    for e in ext:
        w = pywcs.WCS(f[e].header, fobj=f)
        hwcs = w.to_header()
        wcsnamekey = 'WCSNAME' + wkey
        f[e].header.update(key=wcsnamekey, value=wcsname)
        if w.wcs.has_cd():
            pc2cd(hwcs)
        for k in hwcs.keys():
            key = k[:7] + wkey
            f[e].header.update(key=key, value=hwcs[k])
        #norient = np.rad2deg(np.arctan2(hwcs['CD1_2'],hwcs['CD2_2']))
        #okey = 'ORIENT%s' % wkey
        #f[e].header.update(key=okey, value=norient) 
    closefobj(fname, f)

def restoreWCS(f, ext, wcskey=" ", wcsname=" ", clobber=False):
    """
    Copy a WCS with key "WCSKEY" to a primary WCS
     
    Reads in a WCS defined with wcskey and saves it as the primary WCS.
    If clobber is False, writes out new files whose names are the original 
    names with an attached 3 character string  _'WCSKEY'_. 
    Otherwise overwrites the files. Goes sequentially through the list of extensions
    The WCS is restored from the 'SCI' extension but the primary WCS of all 
    extensions with the same EXTVER are updated.
    

    Parameters
    ----------
    f:       string or pyfits.HDUList object
             a file name or a file object
    ext:     an int, a tuple, a python list of integers or a python list 
             of tuples (e.g.('sci',1)) 
             fits extensions to work with
    wcskey:  a charater
             "A"-"Z" - Used for one of 26 alternate WCS definitions.
             or " " - find a key from WCSNAMe value
    wcsname: string (optional)
             if given and wcskey is " ", will try to restore by WCSNAME value
    clobber: boolean
             A flag to define if the original files should be overwritten
             
    See Also
    --------
    wcsutil.archiveWCS - copy the primary WCS as an alternate WCS
    
    """
    if isinstance(f, str):
        if clobber:
            fobj = pyfits.open(f, mode='update')
        else:
            fobj = pyfits.open(f)
    else: 
        fobj = f
        
    if not parpasscheck(fobj, ext, wcskey, wcsname):
        closefobj(f, fobj)
        return
    
    if isinstance(ext, int) or isinstance(ext, tuple):
        ext = [ext]
    
    if not clobber:
        name = (fobj.filename().split('.fits')[0] + '_%s_' + '.fits') %wcskey
    else:
        name = fobj.filename()

    if wcskey == " ":
        if wcsname.strip():
            wkey = getKeyFromName(fobj[1].header, wcsname)
            if not wkey:
                print 'Could not get a key from wcsname %s .' % wcsname
                closefobj(f, fobj)
                return
    else:
        if wcskey not in wcskeys(fobj[1].header):
            print "Could not find alternate WCS with key %s in this file" % wcskey
            closefobj(f, fobj)
            return
        wkey = wcskey    

    for e in ext:
        try:
            extname = fobj[e].header['EXTNAME'].lower()
        except KeyError:
            continue
        #Restore always from a 'SCI' extension but write it out to 'ERR' and 'DQ'
        if extname == 'sci':
            sciver = fobj[e].header['extver']
            try:
                nwcs = pywcs.WCS(fobj[e].header, fobj=fobj, key=wkey)
            except KeyError:
                print 'restoreWCS: Could not read WCS with key %s in file %s,  \
                extension %d' % (wkey, fobj.filename(), e)
                closefobj(f, fobj)
                return #raise
            hwcs = nwcs.to_header()
            
            if nwcs.wcs.has_cd():
                pc2cd(hwcs, key=wkey)
            for k in hwcs.keys():
                key = k[:-1]
                if key in fobj[e].header.keys():
                    fobj[e].header.update(key=key, value = hwcs[k])
                else:
                    continue
            if wcskey == 'O' and fobj[e].header.has_key('TDDALPHA'):
                fobj[e].header['TDDALPHA'] = 0.0
                fobj[e].header['TDDBETA'] = 0.0
            if fobj[e].header.has_key('ORIENTAT'):
                norient = np.rad2deg(np.arctan2(hwcs['CD1_2'+'%s' %wkey],hwcs['CD2_2'+'%s' %wkey]))
                fobj[e].header.update(key='ORIENTAT', value=norient)
        elif extname in ['err', 'dq', 'sdq', 'time', 'samp']:
            cextver = fobj[e].header['extver']
            if cextver == sciver:
                for k in hwcs.keys():
                    key = k[:-1]
                    fobj[e].header.update(key=key, value = hwcs[k])
                if fobj[e].header.has_key('ORIENTAT'):
                    norient = np.rad2deg(np.arctan2(hwcs['CD1_2'+'%s' %wkey],hwcs['CD2_2'+'%s' %wkey]))
                    fobj[e].header.update(key='ORIENTAT', value=norient)
        else:
            continue
            
    if not clobber:
        fobj.writeto(name)
    closefobj(f, fobj)

def deleteWCS(fname, ext, wcskey=" ", wcsname=" "):
    """
    Delete an alternate WCS defined with wcskey.
    If wcskey is " " try to get a key from WCSNAME.
    
    Parameters
    ----------
    fname:   sting or a pyfits.HDUList object
    ext:     an int, a tuple, a python list of integers or a python list of tuples (e.g.('sci',1))
             fits extensions to work with
    wcskey:  one of 'A'-'Z' or " "
    wcsname: string
             Name of alternate WCS description
    """
    if isinstance(fname, str):
        fobj = pyfits.open(fname, mode='update')
    else:
        fobj = fname
    
    if not parpasscheck(fobj, ext, wcskey, wcsname):
        closefobj(fname, fobj)
        return
    
    if isinstance(ext, int) or isinstance(ext, tuple):
        ext = [ext]

    # Do not allow deleting the original WCS.
    if wcskey == 'O':
        print "Wcskey 'O' is reserved for the original WCS and should not be deleted." 
        closefobj(fname, fobj)
        return
    
    if wcskey == " ": 
        # try getting the key from WCSNAME
        if wcsname == " ":
            print "Could not get a valid key from header"
            closefobj(fname, fobj)
            return
        else:
            wkey = getKeyFromName(fobj[1].header, wcsname)
            if not wkey:
                print 'Could not get a key: wcsname "%s" not found in header.' % wcsname
                closefobj(fname, fobj)
                return
    else:
        if wcskey not in wcskeys(fobj[1].header):
            print "Could not find alternate WCS with key %s in this file" % wcskey
            closefobj(fname, fobj)
            return
        wkey = wcskey    
        
    prexts = []
    for i in ext:
        hdr = fobj[i].header
        try:
            w = pywcs.WCS(hdr, fobj, key=wkey)
        except KeyError:
            continue
        hwcs = w.to_header()
        if w.wcs.has_cd():
            pc2cd(hwcs, key=wkey)
        for k in hwcs.keys():
            del hdr[k]
            #del hdr['ORIENT'+wkey]
        prexts.append(i)
    if prexts != []:
        print 'Deleted all instances of WCS with key %s in extensions' % wkey, prexts
    else:
        print "Did not find WCS with key %s in any of the extensions" % wkey
    closefobj(fname, fobj)

#header operations
def wcskeys(header):
    """
    Returns a list of characters used in the header for alternate 
    WCS description with WCSNAME keyword
    
    Parameters
    ----------
    hdr: pyfits.Header
    """
    assert isinstance(header, pyfits.Header), "Requires a pyfits.Header object as input"
    names = header["WCSNAME*"]
    return [key.split('WCSNAME')[1].upper() for key in names.keys()]

def wcsnames(header):
    """
    Returns a dictionary of wcskey: WCSNAME pairs
    
    Parameters
    ----------
    header: pyfits.Header
    """
    names = header["WCSNAME*"]
    d = {}
    for c in names:
        d[c.key[-1]] = c.value
    return d
    
def available_wcskeys(header):
    """
    Returns a list of characters which are not used in the header
    with WCSNAME keyword. Any of them can be used to save a new
    WCS.
    
    Parameters
    ----------
    header: pyfits.Header
    """
    assert isinstance(header, pyfits.Header), "Requires a pyfits.Header object as input"
    all_keys = list(string.ascii_uppercase)
    used_keys = wcskeys(header)
    try: 
        used_keys.remove("")
    except ValueError: 
        pass
    [all_keys.remove(key) for key in used_keys]
    return all_keys

def next_wcskey(header):
    """
    Returns next available character to be used for an alternate WCS
    
    Parameters
    ----------
    header: pyfits.Header    
    """
    assert isinstance(header, pyfits.Header), "Requires a pyfits.Header object as input"
    allkeys = available_wcskeys(header)
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

def parpasscheck(fobj, ext, wcskey, wcsname, clobber=True):

    if not isinstance(fobj,pyfits.HDUList):
        print "First parameter must be a fits file object or a file name."
        return False
    try:
        assert (fobj.fileinfo(0)['filemode'] == 'update')
    except AssertionError:
        print "First parameter must be a file name or a file object opened in 'update' mode."
        return False

    if not isinstance(ext, int) and not isinstance(ext, tuple) \
       and not isinstance(ext, list):
        print "Ext must be integer, tuple, a list of int extension numbers, \
        or a list of tuples representing a fits extension, for example ('sci', 1)."
        return False
    
    if len(wcskey) != 1:
        print 'Parameter wcskey must be a character - one of "A"-"Z" or " "'
        return False
    
    if wcskey == " ": 
        # try getting the key from WCSNAME
        """
        if wcsname == " " or wcsname == "":
            #wkey = next_wcskey(f[1].header)
            #if not wkey:
            #    print "Could not get a valid key from header"
            return False
        """
        if wcsname.strip():
            wkey = getKeyFromName(fobj[1].header, wcsname)
            if wkey and not clobber:
                print 'Wcsname %s is already used.' % wcsname
                print 'Use "wcsutil.next_wcskey" to obtain a valid wcskey'
                print 'or use "clobber=True" to overwrite the values.'
                return False
    else:
        if wcskey not in available_wcskeys(fobj[1].header):
            if clobber==False:
                print 'Wcskey %s is already used.' % wcskey
                print 'Use "wcsutil.next_wcskey" to obtain a valid wcskey'
                print 'or use "clobber=True" to overwrite the values.'
                return False


    return True

def closefobj(fname,f):
    """
    Functions in this module accept as input a file name or a file object.
    If the input was a file name (string) we close the object. If the user
    passed a file object we leave it to the user to close it.
    """
    if isinstance(fname, str):
        f.close()
        

    