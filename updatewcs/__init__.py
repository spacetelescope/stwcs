import os
import pyfits
#from .. wcsutil import HSTWCS
from hstwcs.wcsutil import HSTWCS
from hstwcs.mappings import allowed_corrections
#from .. mappings import allowed_corrections
import corrections, makewcs
import dgeo
import time
from pytools import parseinput, fileutil

#NB! the order of corrections matters

__docformat__ = 'restructuredtext'


def updatewcs(input, vacorr=True, tddcorr=True, checkfiles=True):
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
              If True, time dependen correction will be applied to the distortion model
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
        instr = pyfits.getval(f, 'INSTRUME')
        try:
            acorr = setCorrections(instr,vacorr=vacorr, tddcorr=tddcorr)
        except KeyError:
            print 'Unsupported instrument %s ' %instr
            print 'Removing %s from list of processed files\n' % f
            files.remove(f)
            continue
        
        makecorr(f, acorr)
    return files

def makecorr(fname, acorr):
    """
    Purpose
    =======
    Applies corrections to a single file
    
    :Parameters:
    `fname`: string
             file name
    `acorr`: list
             list of corrections to be applied
            
    """
    f = pyfits.open(fname, mode='update')
    nrefchip, nrefext = getNrefchip(f)
    primhdr = f[0].header
    
    
    for extn in f:
        refwcs = HSTWCS(primhdr, f[nrefext].header)
        refwcs.archive_kw()
        refwcs.readModel()
        if extn.header.has_key('extname') and extn.header['extname'].lower() == 'sci':
            hdr = extn.header
            owcs = HSTWCS(primhdr, hdr)
            owcs.archive_kw()
            owcs.readModel()
            for c in acorr:
                owcs.__setattr__('DO'+c, 'PERFORM')
                corr_klass = corrections.__getattribute__(c)
                corr_klass(owcs, refwcs)
            
            
    
    #always do dgeo correction
    if applyDgeoCorr(fname):
        dgeo.DGEO(f)
    f.close()

        
def setCorrections(instrument, vacorr=True, tddcorr=True):
    """
    Purpose
    =======
    Creates a list of corrections to be applied to a file.
    based on user input paramters and allowed corrections
    for the instrument, which are defined in mappings.py.
    """
    acorr = allowed_corrections[instrument]
    if 'VACorr' in acorr and not vacorr:  acorr.remove('VACorr')
    if 'TDDCorr' in acorr and not tddcorr: acorr.remove('TDDCorr')
    if 'DGEOCorr' in acorr and not dgeocorr: acorr.remove('DGEOCorr')

    return acorr



def applyDgeoCorr(fname):
    """
    Purpose
    =======
    Adds dgeo extensions to files based on the DGEOFILE keyword in the primary 
    header. This is a default correction and will always run in the pipeline.
    The file used to generate the extensions is 
    recorded in the DGEOFILE keyword in each science extension.
    If 'DGEOFILE' in the primary header is different from 'DGEOFILE' in the 
    extension header and the file exists on disk and is a 'new type' dgeofile, 
    then the dgeo extensions will be updated.
    """
    applyDGEOCorr = True
    try:
        # get DGEOFILE kw from primary header
        fdgeo0 = pyfits.getval(fname, 'DGEOFILE')
        fdgeo0 = fileutil.osfn(fdgeo0)
        if not fileutil.findFile(fdgeo0):
            print 'Kw DGEOFILE exists in primary header but file %s not found\n' % fdgeo0
            print 'DGEO correction will not be applied\n'
            applyDGEOCorr = False
            return applyDGEOCorr 
        try:
            # get DGEOFILE kw from first extension header
            fdgeo1 = pyfits.getval(fname, 'DGEOFILE', ext=1)
            fdgeo1 = fileutil.osfn(fdgeo1)
            if fdgeo1 and fileutil.findFile(fdgeo1):
                if fdgeo0 != fdgeo1:
                    applyDGEOCorr = True
                else:
                    applyDGEOCorr = False
            else: 
                # dgeo file defined in first extension may not be found
                # but if a valid kw exists in the primary header, dgeo should be applied.
                applyDGEOCorr = True
        except KeyError:
            # the case of DGEOFILE kw present in primary header but missing 
            # in first extension header
            applyDGEOCorr = True
    except KeyError:

        print 'DGEOFILE keyword not found in primary header'
        applyDGEOCorr = False
    
    if isOldStyleDGEO(fname, fdgeo0):
        applyDGEOCorr = False
        
    return applyDGEOCorr

def isOldStyleDGEO(fname, dgname):
    # checks if the file defined in a DGEOFILE kw is a full size 
    # (old style) image
    
    sci_naxis1 = pyfits.getval(fname, 'NAXIS1', ext=1)
    sci_naxis2 = pyfits.getval(fname, 'NAXIS2', ext=1)
    dg_naxis1 = pyfits.getval(dgname, 'NAXIS1', ext=1)
    dg_naxis2 = pyfits.getval(dgname, 'NAXIS2', ext=1)
    if sci_naxis1 <= dg_naxis1 or sci_naxis2 <= dg_naxis2:
        print 'Only full size (old style) XY file was found.'
        print 'DGEO correction will not be applied.\n'
        return True
    else:
        return False
    
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
    return Nrefchip, Nrefext

def checkFiles(input):
    """
    Purpose
    =======
    Checks that input files are in the correct format.
    Converts geis and waiver fits files to multietension fits.
    """
    from pytools.check_files import geis2mef, waiver2mef
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
    return newfiles

