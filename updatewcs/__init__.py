from __future__ import division # confidence high

import os
import pyfits
from stwcs.wcsutil import HSTWCS

from stwcs import utils
import corrections, makewcs
import dgeo, det2im
from pytools import parseinput, fileutil
import apply_corrections

#Note: The order of corrections is important

__docformat__ = 'restructuredtext'

__version__ = '0.5'

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
        
        #restore the original WCS keywords
        utils.restoreWCS(f, wcskey='O', clobber=True)
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
    #Determine the reference chip and create the reference HSTWCS object
    nrefchip, nrefext = getNrefchip(f)
    ref_wcs = HSTWCS(fobj=f, ext=nrefext)
    ref_wcs.readModel(update=True,header=f[nrefext].header)
    ref_wcs.copyWCS(header=f[nrefext].header, wcskey='O', wcsname='OPUS', clobber=True)
    
    if 'DET2IMCorr' in allowed_corr:
        det2im.DET2IMCorr.updateWCS(f)
        
    for i in range(len(f))[1:]:
        # Perhaps all ext headers should be corrected (to be consistent)
        extn = f[i]
        
        if extn.header.has_key('extname'):
            extname = extn.header['extname'].lower()
            if  extname == 'sci':
                
                sciextver = extn.header['extver']
                ref_wcs.restore(f[nrefext].header, wcskey="O")
    
                hdr = extn.header
                ext_wcs = HSTWCS(fobj=f, ext=i)
                ext_wcs.copyWCS(header=hdr, wcskey='O', wcsname='OPUS', clobber=True)
                ext_wcs.readModel(update=True,header=hdr)
                for c in allowed_corr:
                    if c != 'DGEOCorr' and c != 'DET2IMCorr':
                        corr_klass = corrections.__getattribute__(c)
                        kw2update = corr_klass.updateWCS(ext_wcs, ref_wcs)
                        for kw in kw2update:
                            hdr.update(kw, kw2update[kw])
                    
                if wkey is not None:
                    # archive the updated primary WCS
                    if wkey == " " :
                        idcname = os.path.split(fileutil.osfn(ext_wcs.idctab))[1]
                        wname = ''.join(['IDC_',idcname.split('_idc.fits')[0]])
                        wkey = getKey(hdr, wname)
                        #in this case clobber = true, to allow updatewcs to be run repeatedly
                        ext_wcs.copyWCS(header=hdr, wcskey=wkey, wcsname=wname, clobber=True)
                    else:
                        #clobber is set to False as a warning to users
                        ext_wcs.copyWCS(header=hdr, wcskey=wkey, wcsname=wname, clobber=False)
                    
            elif extname in ['err', 'dq', 'sdq']:
                cextver = extn.header['extver']
                if cextver == sciextver:
                    ext_wcs.copyWCS(header=extn.header, wcskey=" ", wcsname=" ")
            else:
                cextver = extn.header['extver']
                continue
    
    if 'DGEOCorr' in allowed_corr:
        kw2update = dgeo.DGEOCorr.updateWCS(f)
        for kw in kw2update:
            f[1].header.update(kw, kw2update[kw])       
            
    f.close()

def getKey(header, wcsname):
    """
    If WCSNAME is found in header, return its key, else return 
    the next available key. This is used to update a specific WCS
    repeatedly and not generate new keys every time.
    """
    wkey = utils.next_wcskey(header)
    names = utils.wcsnames(header)
    for item in names.items():
        if item[1] == wcsname:
            wkey = item[0]
    return wkey

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


def getCorrections(instrument):
    """
    Print corrections available for an instrument
    
    :Parameters:
    `instrument`: string, one of 'WFPC2', 'NICMOS', 'STIS', 'ACS', 'WFC3'
    """
    acorr = apply_corrections.allowed_corrections[instrument]
    
    print "The following corrections will be performed for instrument %s\n" % instrument
    for c in acorr: print c,': ' ,  apply_corrections.cnames[c]
    