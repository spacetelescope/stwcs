import os
import pyfits
#from .. wcsutil import HSTWCS
from updatewcs.wcsutil import HSTWCS

#from .. mappings import allowed_corrections
from updatewcs import utils
import corrections, makewcs
import dgeo
import time
from pytools import parseinput, fileutil
import apply_corrections

#Note: The order of corrections is important

__docformat__ = 'restructuredtext'

__version__ = '0.3'

def updatewcs(input, vacorr=True, tddcorr=True, dgeocorr=True, checkfiles=True):
    """
    Purpose
    =======
    Applies corrections to the WCS keywords.
    
    Example
    =======
    >>>from hstwcs import updatewcs
    >>>updatewcs.updatewcs(filename)
    
    Dependencies 
    ============
    `pytools`
    `pyfits`
    `pywcs`
    `numpy`

    :Parameters:
    `input`: a python list of file names or a string (wild card characters allowed)
             input files may be in fits, geis or waiver fits format
    `vacorr`: boolean
              If True, vecocity aberration correction will be applied
    `tddcorr`: boolean
              If True, time dependent distortion correction will be applied 
    `checkfiles`: boolean
              If True, the format of the input files will be checked,
              geis and waiver fits files will be converted to MEF format.
              Default value is True for standalone mode.
    """
    
    files = parseinput.parseinput(input)[0]
    if checkfiles:
        files = checkFiles(files)
        if not files:
            print 'No valid input, quitting ...\n'
            return
    for f in files:
        acorr = apply_corrections.setCorrections(f, vacorr=vacorr, tddcorr=tddcorr,dgeocorr=dgeocorr)
        #restore the original WCS keywords
        utils.restoreWCS(f)
        makecorr(f, acorr)
    return files

def makecorr(fname, allowed_corr):
    """
    Purpose
    =======
    Applies corrections to the WCS of a single file
    
    :Parameters:
    `fname`: string
             file name
    `acorr`: list
             list of corrections to be applied
            
    """
    f = pyfits.open(fname, mode='update')
    #Determine the reference chip and create the reference HSTWCS object
    nrefchip, nrefext = getNrefchip(f)
    ref_wcs = HSTWCS(fobj=f, ext=nrefext)
    ref_wcs.readModel(update=True,header=f[nrefext].header)
    utils.write_archive(f[nrefext].header)
            
    for i in range(len(f))[1:]:
        # Perhaps all ext headers should be corrected (to be consistent)
        extn = f[i]
        if extn.header.has_key('extname') and extn.header['extname'].lower() == 'sci':
            ref_wcs.restore(f[nrefext].header)
            hdr = extn.header
            utils.write_archive(hdr)
            ext_wcs = HSTWCS(fobj=f, ext=i)
            ext_wcs.readModel(update=True,header=hdr)
            for c in allowed_corr:
                if c != 'DGEOCorr':
                    corr_klass = corrections.__getattribute__(c)
                    kw2update = corr_klass.updateWCS(ext_wcs, ref_wcs)
                    for kw in kw2update:
                        hdr.update(kw, kw2update[kw])
        
    if 'DGEOCorr' in allowed_corr:
        kw2update = dgeo.DGEOCorr.updateWCS(f)
        for kw in kw2update:
            f[1].header.update(kw, kw2update[kw])
            
    f.close()
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
        if detector == 'WCS':
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

