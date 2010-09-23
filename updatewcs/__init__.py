from __future__ import division # confidence high

import os
import pyfits
import numpy as np
from stwcs import wcsutil
from stwcs.wcsutil import HSTWCS
import pywcs

import utils, corrections, makewcs
import dgeo, det2im
from pytools import parseinput, fileutil
import apply_corrections

#Note: The order of corrections is important

__docformat__ = 'restructuredtext'

__version__ = '0.8'

def updatewcs(input, vacorr=True, tddcorr=True, dgeocorr=True, d2imcorr=True, 
              checkfiles=True, wcskey=" ", wcsname=" ", clobber=False):
    """
    Purpose
    =======
    Updates HST science files with the best available calibration information.
    This allows users to retrieve from the archive self contained science files 
    which do not require additional reference files.
    
    Basic WCS keywords are updated in the process and new keywords (following WCS 
    Paper IV and the SIP convention) as well as new extensions are added to the science files.
    
    
    Example
    =======
    >>>from stwcs import updatewcs
    >>>updatewcs.updatewcs(filename)
    
    Dependencies 
    ============
    `pytools`
    `pyfits`
    `pywcs`

    :Parameters:
    `input`: a python list of file names or a string (wild card characters allowed)
             input files may be in fits, geis or waiver fits format
    `vacorr`: boolean
              If True, vecocity aberration correction will be applied
    `tddcorr`: boolean
              If True, time dependent distortion correction will be applied 
    `dgeocorr`: boolean
              If True, a Lookup table distortion will be applied
    `d2imcorr`: boolean
              If True, detector to image correction will be applied
    `checkfiles`: boolean
              If True, the format of the input files will be checked,
              geis and waiver fits files will be converted to MEF format.
              Default value is True for standalone mode.
    `wcskey`: None, one character string A-Z or an empty string of length 1
              If None - the primary WCS is not archived
              If an empty string - the next available wcskey is used for the archive
              A-Z - use this key to archive the WCS
    `wcsname`: a string
              The name under which the primary WCS is archived after it is updated.
              If an empty string (default), the name of the idctable is used as 
              a base.
    `clobber`: boolean
              a flag for reusing the wcskey when archiving the primary WCS
    """
    files = parseinput.parseinput(input)[0]
    if checkfiles:
        files = checkFiles(files)
        if not files:
            print 'No valid input, quitting ...\n'
            return 
    for f in files:
        acorr = apply_corrections.setCorrections(f, vacorr=vacorr, \
            tddcorr=tddcorr,dgeocorr=dgeocorr, d2imcorr=d2imcorr)
        
        if 'MakeWCS' in acorr and newIDCTAB(f):
            print "New IDCTAB file detected. This invalidates all WCS's." 
            print "Deleting all previous WCS's"
            cleanWCS(f)
            
        #restore the original WCS keywords
        #wcsutil.restoreWCS(f, ext=[], wcskey='O', clobber=True)
        makecorr(f, acorr, wkey=wcskey, wname=wcsname, clobber=False)
    return files

def makecorr(fname, allowed_corr, wkey=" ", wname=" ", clobber=False):
    """
    Purpose
    =======
    Applies corrections to the WCS of a single file
    
    :Parameters:
    `fname`: string
             file name
    `acorr`: list
             list of corrections to be applied
    `wkey`: None, one character string A-Z or an empty string of length 1
              If None - the primary WCS is not archived
              If an empty string - the next available wcskey is used for the archive
              A-Z - use this key to archive the WCS
    `wname`: a string
              The name under which the primary WCS is archived after it is updated.
              If an empty string (default), the name of the idctable is used as 
              a base.
    `clobber`: boolean
              a flag for reusing the wcskey when archiving the primary WCS
    """
    f = pyfits.open(fname, mode='update')
    #restore the original WCS keywords
    #wcsutil.restoreWCS(f, ext=[], wcskey='O', clobber=True)
    #Determine the reference chip and create the reference HSTWCS object
    nrefchip, nrefext = getNrefchip(f)
    wcsutil.restoreWCS(f, nrefext, wcskey='O', clobber=True)
    rwcs = HSTWCS(fobj=f, ext=nrefext)
    rwcs.readModel(update=True,header=f[nrefext].header)
    
    wcsutil.archiveWCS(f, nrefext, 'O', wcsname='OPUS', clobber=True)
    
    if 'DET2IMCorr' in allowed_corr:
        det2im.DET2IMCorr.updateWCS(f)
    
    # get a wcskey and wcsname from the first extension header
    idcname = fileutil.osfn(rwcs.idctab)
    key, name = getKeyName(f[1].header, wkey, wname, idcname)
    
    for i in range(len(f))[1:]:
        extn = f[i]
        
        if extn.header.has_key('extname'):
            extname = extn.header['extname'].lower()
            if  extname == 'sci':
                wcsutil.restoreWCS(f, ext=i, wcskey='O', clobber=True)
                sciextver = extn.header['extver']
                ref_wcs = rwcs.deepcopy()
                hdr = extn.header
                ext_wcs = HSTWCS(fobj=f, ext=i)
                wcsutil.archiveWCS(f, ext=i, wcskey="O", wcsname="OPUS", clobber=True)
                ext_wcs.readModel(update=True,header=hdr)
                for c in allowed_corr:
                    if c != 'DGEOCorr' and c != 'DET2IMCorr':
                        corr_klass = corrections.__getattribute__(c)
                        kw2update = corr_klass.updateWCS(ext_wcs, ref_wcs)
                        for kw in kw2update:
                            hdr.update(kw, kw2update[kw])
                #if wkey is None, do not archive the primary WCS  
                if key is not None:
                    wcsutil.archiveWCS(f, ext=i, wcskey=key, wcsname=name, clobber=True)
            elif extname in ['err', 'dq', 'sdq', 'samp', 'time']:
                cextver = extn.header['extver']
                if cextver == sciextver:
                    hdr = f[('SCI',sciextver)].header
                    w = pywcs.WCS(hdr, f)
                    copyWCS(w, extn.header, key, name)
            else:
                continue
    
    if 'DGEOCorr' in allowed_corr:
        kw2update = dgeo.DGEOCorr.updateWCS(f)
        for kw in kw2update:
            f[1].header.update(kw, kw2update[kw])       
        
    f.close()

def getKeyName(hdr, wkey, wname, idcname):
    if wkey is not None: # archive the primary WCS
        if wkey == " ":
            if wname == " " :
                # get the next available key and use the IDCTABLE name as WCSNAME
                idcname = os.path.split(idcname)[1]
                name = ''.join(['IDC_',idcname.split('_idc.fits')[0]])
                key = wcsutil.getKeyFromName(hdr, name)
                if not key:
                    key = wcsutil.next_wcskey(hdr)
            else:
                #try to get a key from WCSNAME
                # if not - get the next availabble key
                name = wname
                key = wcsutil.getKeyFromName(hdr, wname)
                if not wkey:
                    key = wcsutil.next_wcskey(hdr)
        else:
            key = wkey
            name = wname
    return key, name

def copyWCS(w, hdr, wkey, wname):
    """
    This is a convenience function to copy a WCS object 
    to a header as a primary WCS. It is used only to copy the 
    WCS of the 'SCI' extension to the headers of 'ERR', 'DQ', 'SDQ',
    'TIME' or 'SAMP' extensions.
    """
    hwcs = w.to_header()
    
    if w.wcs.has_cd():
        wcsutil.pc2cd(hwcs)
    for k in hwcs.keys():
        key = k+wkey
        hdr.update(key=key, value=hwcs[k])
    norient = np.rad2deg(np.arctan2(hwcs['CD1_2'],hwcs['CD2_2']))
    okey = 'ORIENT%s' % wkey
    hdr.update(key=okey, value=norient) 
    
def getNrefchip(fobj):
    """
    This handles the fact that WFPC2 subarray observations
    may not include chip 3 which is the default reference chip for
    full observations. Also for subarrays chip 3  may not be the third
    extension in a MEF file. 
    """
    Nrefext = 1
    instrument = fobj[0].header['INSTRUME']
    if instrument == 'WFPC2':
        detectors = [img.header['DETECTOR'] for img in fobj[1:]]

        if 3 not in detectors:
            Nrefchip=detectors[0]
            Nrefext = 1
        else:
            Nrefchip = 3
            Nrefext = detectors.index(3) + 1
    elif instrument == 'ACS':
        detector = fobj[0].header['DETECTOR']
        if detector == 'WFC':
            Nrefchip =2
        else:
            Nrefchip = 1
    elif instrument == 'NICMOS':
        Nrefchip = fobj[0].header['CAMERA']
    elif instrument == 'WFC3':
        detector = fobj[0].header['DETECTOR']
        if detector == 'UVIS':
            Nrefchip =2
        else:
            Nrefchip = 1
    else:
        Nrefchip = 1
    return Nrefchip, Nrefext

def checkFiles(input):
    """
    Purpose
    =======
    Checks that input files are in the correct format.
    Converts geis and waiver fits files to multiextension fits.
    """
    from pytools.check_files import geis2mef, waiver2mef, checkFiles
    removed_files = []
    newfiles = []
    for file in input:
        try:
                imgfits,imgtype = fileutil.isFits(file)
        except IOError:
            print "Warning:  File %s could not be found\n" %file
            print "Removing file %s from input list" %file
            removed_files.append(file)
            continue
        # Check for existence of waiver FITS input, and quit if found.
        # Or should we print a warning and continue but not use that file
        if imgfits: 
            if imgtype == 'waiver':
                newfilename = waiver2mef(file, convert_dq=True)
                if newfilename == None:
                    print "Removing file %s from input list - could not convert waiver to mef" %file
                    removed_files.append(file)
                else:
                    newfiles.append(newfilename)
            else:
                newfiles.append(file)

        # If a GEIS image is provided as input, create a new MEF file with 
        # a name generated using 'buildFITSName()'
        # Convert the corresponding data quality file if present    
        if not imgfits:
            newfilename = geis2mef(file, convert_dq=True)
            if newfilename == None:
                print "Removing file %s from input list - could not convert geis to mef" %file
                removed_files.append(file)
            else:
                newfiles.append(newfilename)
    if removed_files:
        print 'The following files will be removed from the list of files to be processed :\n'
        for f in removed_files:
            print f
            
    newfiles = checkFiles(newfiles)[0]
    
    return newfiles

def newIDCTAB(fname):
    #When this is called we know there's a kw IDCTAB in the header
    idctab = fileutil.osfn(pyfits.getval(fname, 'IDCTAB'))
    try:
        #check for the presence of IDCTAB in the first extension
        oldidctab = fileutil.osfn(pyfits.getval(fname, 'IDCTAB', ext=1))
    except KeyError:
        return False
    if idctab == oldidctab:
        return False
    else:
        return True
    
def cleanWCS(fname):
    # A new IDCTAB means all previously computed WCS's are invalid
    # We are deleting all of them except the original OPUS WCS.nvalidates all WCS's.
    keys = wcsutil.wcskeys(pyfits.getheader(fname, ext=1))
    f = pyfits.open(fname, mode='update')
    fext = range(len(f))
    for key in keys:
        wcsutil.deleteWCS(fname, ext=fext,wcskey=key)
            
def getCorrections(instrument):
    """
    Print corrections available for an instrument
    
    :Parameters:
    `instrument`: string, one of 'WFPC2', 'NICMOS', 'STIS', 'ACS', 'WFC3'
    """
    acorr = apply_corrections.allowed_corrections[instrument]
    
    print "The following corrections will be performed for instrument %s\n" % instrument
    for c in acorr: print c,': ' ,  apply_corrections.cnames[c]
    