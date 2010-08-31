from __future__ import division # confidence high

from pytools import parseinput, fileutil
import pyfits
import pywcs
import numpy as np
import string 

def restoreWCS(fnames, wcskey, clobber=False):
    """
    Purpose
    =======
    Reads in a WCS defined with wcskey and saves it as the primary WCS.
    If clobber is False, writes out new files whose names are the original 
    names with an attached 3 character string representing _'WCSKEY'_. 
    Otherwise overwrites the files. The WCS is restored from the 'SCI' 
    extension but the primary WCS of all extensions are updated.
    
    WCS keywords:
    'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 
    'CRVAL*',
    'CTYPE*',
    'CRPIX*', 
    'CDELT*',
    'CUNIT*',
    'ORIENTAT' - ?
    'TDDALPHA', 'TDDBETA'
    'A_x_x', B_x_x' - SIP coefficients
    'CPERROR*', 'CPDIS*', 'DP*', 
    
    `fnames`: a python list of file names, a string of comma separated file names, 
              an @file
    `wcskey`: a charater
              Used for one of 26 alternate WCS definitions.
    `clobber`: boolean
              A flag to define if the original files should be overwritten
    """
    files = parseinput.parseinput(fnames)[0]
    for f in files:
        isfits, ftype = fileutil.isFits(f)
        if not isfits or (isfits and ftype == 'waiver'):
            print "RestoreWCS works only with true fits files."
            return
        else:
            if clobber:
                print 'Overwriting original files\n'
                fobj = pyfits.open(f, mode='update')
                name = f
            else:
                fobj = pyfits.open(f)
                name = (f.split('.fits')[0] + '_%s_' + '.fits') %wcskey
            for e in range(len(fobj)):
                try:
                    extname = fobj[e].header['EXTNAME'].lower()
                except KeyError:
                    continue
                #Restore always from a 'SCI' extension but write it out to 'ERR' and 'DQ'
                if extname == 'sci':
                    sciver = fobj[e].header['extver']
                    try:
                        nwcs = pywcs.WCS(fobj[e].header, fobj=fobj, key=wcskey)
                    except:
                        print 'utils.restoreWCS: Could not read WCS with key %s in file %s,  \
                        extension %d' % (wcskey, f, e)
                        return #raise
                    hwcs = nwcs.to_header()
                    if nwcs.wcs.has_pc():
                        for c in ['1_1', '1_2', '2_1', '2_2']:
                            del hwcs['CD'+c+wcskey]
                    elif nwcs.wcs.has_cd():
                        for c in ['1_1', '1_2', '2_1', '2_2']:
                            hwcs.update(key='CD'+c+wcskey, value=hwcs['PC'+c+wcskey])
                            del hwcs['PC'+c]
                    for k in hwcs.keys():
                        key = k[:-1]
                        if key in fobj[e].header.keys():
                            fobj[e].header.update(key=key, value = hwcs[k])
                        else:
                            continue
                    if wcskey == 'O':
                        fobj[e].header['TDDALPHA'] = 0.0
                        fobj[e].header['TDDBETA'] = 0.0
                    if fobj[e].header.has_key('ORIENTAT'):
                        cd12 = 'CD1_2%s' % wcskey
                        cd22 = 'CD2_2%s' % wcskey
                        norient = np.rad2deg(np.arctan2(hwcs[cd12],hwcs[cd22]))
                        fobj[e].header.update(key='ORIENTAT', value=norient)
                elif extname in ['err', 'dq', 'sdq']:
                    cextver = fobj[e].header['extver']
                    if cextver == sciver:
                        for k in hwcs.keys():
                            key = k[:-1]
                            fobj[e].header.update(key=key, value = hwcs[k])
                        if fobj[e].header.has_key('ORIENTAT'):
                            cd12 = 'CD1_2%s' % wcskey
                            cd22 = 'CD2_2%s' % wcskey
                            norient = np.rad2deg(np.arctan2(hwcs[cd12],hwcs[cd22]))
                            fobj[e].header.update(key='ORIENTAT', value=norient)
                else:
                    continue
                
            if not clobber:
                fobj.writeto(name)
            fobj.close()

def archiveWCS(fname, ext, wcskey, wcsname=" "):
    """
    Copy the primary WCS to an alternate WCS 
    with wcskey and name WCSNAME.
    """
    f = pyfits.open(fname, mode='update')
    w = pywcs.WCS(f[ext].header, fobj=f)
    assert len(wcskey) == 1
    if wcskey == " ": 
        print "Please provide a valid wcskey for this WCS."
        print 'Use "utils.next_wcskey" to obtain a valid wcskey.'
        print 'Use utils.restoreWCS to write to the primary WCS.'
        return
    if wcskey not in available_wcskeys(f[ext].header):
        print 'wcskey %s is already used in this header.' % wcskey
        print 'Use "utils.next_wcskey" to obtain a valid wcskey'
    hwcs = w.to_header()
    wkey = 'WCSNAME' + wcskey
    f[ext].header.update(key=wkey, value=wcsname)
    if w.wcs.has_pc():
        for c in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
            del hwcs[c]
    elif w.wcs.has_cd():
        for c in ['PC1_1', 'PC1_2', 'PC2_1', 'PC2_2']:
            del hwcs[c]
    for k in hwcs.keys():
        key = k+wcskey
        f[ext].header.update(key=key, value = hwcs[k])
    norient = np.rad2deg(np.arctan2(hwcs['CD1_2'],hwcs['CD2_2']))
    okey = 'ORIENT%s' % wcskey
    f[ext].header.update(key=okey, value=norient) 
    f.close()
    

def diff_angles(a,b):
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
        binned = fobj['SCI', extver].header.get('BINAXIS',1)
    return binned

def wcsnames(header):
    """
    Purpose
    =======
    Return a dictionary of wcskey: WCSNAME pairs
    """
    names = header["WCSNAME*"]
    d = {}
    for c in names:
        d[c.key[-1]] = c.value
    return d
    

def wcskeys(header):
    """
    Purpose
    =======
    Returns a list of characters used in the header for alternate 
    WCS description via WCSNAME keyword
    
    `header`: pyfits.Header
    """
    names = header["WCSNAME*"]
    return [key.split('WCSNAME')[1] for key in names.keys()]

def available_wcskeys(header):
    """
    Purpose
    =======
    Returns a list of characters which are not used in the header
    with WCSNAME keyword. Any of them can be used to save a new
    WCS.
    
    `header`: pyfits.Header
    """
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
    Purpose
    =======
    Returns next available character to be used for an alternate WCS
    
    `header`: pyfits.Header    
    """
    allkeys = available_wcskeys(header)
    if allkeys != []:
        return allkeys[0]
    else:
        return None
    