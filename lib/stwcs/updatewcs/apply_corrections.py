from __future__ import division # confidence high

import os
import pyfits
import time
from stsci.tools import fileutil
import os.path
from stwcs.wcsutil import altwcs

import logging
logger = logging.getLogger("stwcs.updatewcs.apply_corrections")

#Note: The order of corrections is important

__docformat__ = 'restructuredtext'

# A dictionary which lists the allowed corrections for each instrument.
# These are the default corrections applied also in the pipeline.

allowed_corrections={'WFPC2': ['DET2IMCorr', 'MakeWCS','CompSIP', 'VACorr'],
                    'ACS': ['DET2IMCorr', 'TDDCorr', 'MakeWCS', 'CompSIP','VACorr', 'NPOLCorr'],
                    'STIS': ['MakeWCS', 'CompSIP','VACorr'],
                    'NICMOS': ['MakeWCS', 'CompSIP','VACorr'],
                    'WFC3': ['DET2IMCorr', 'MakeWCS', 'CompSIP','VACorr', 'NPOLCorr'],
                    }

cnames = {'DET2IMCorr': 'Detector to Image Correction',
         'TDDCorr': 'Time Dependent Distortion Correction',
         'MakeWCS': 'Recalculate basic WCS keywords based on the distortion model',
         'CompSIP': 'Given IDCTAB distortion model calculate the SIP coefficients',
         'VACorr':  'Velocity Aberration Correction',
         'NPOLCorr': 'Lookup Table Distortion'
         }

def setCorrections(fname, vacorr=True, tddcorr=True, npolcorr=True, d2imcorr=True):
    """
    Creates a list of corrections to be applied to a file
    based on user input paramters and allowed corrections
    for the instrument.
    """
    instrument = pyfits.getval(fname, 'INSTRUME')
    # make a copy of this list !
    acorr = allowed_corrections[instrument][:]

    # Check if idctab is present on disk
    # If kw IDCTAB is present in the header but the file is
    # not found on disk, do not run TDDCorr, MakeCWS and CompSIP
    if not foundIDCTAB(fname):
        if 'TDDCorr' in acorr: acorr.remove('TDDCorr')
        if 'MakeWCS' in acorr: acorr.remove('MakeWCS')
        if 'CompSIP' in acorr: acorr.remove('CompSIP')

    if 'VACorr' in acorr and vacorr==False:  acorr.remove('VACorr')
    if 'TDDCorr' in acorr:
        tddcorr = applyTDDCorr(fname, tddcorr)
        if tddcorr == False: acorr.remove('TDDCorr')

    if 'NPOLCorr' in acorr:
        npolcorr = applyNpolCorr(fname, npolcorr)
        if npolcorr == False: acorr.remove('NPOLCorr')
    if 'DET2IMCorr' in acorr:
        d2imcorr = applyD2ImCorr(fname, d2imcorr)
        if d2imcorr == False: acorr.remove('DET2IMCorr')
    logger.info("\n\tCorrections to be applied to %s: %s " % (fname, str(acorr)))
    return acorr

def foundIDCTAB(fname):
    try:
        idctab = fileutil.osfn(pyfits.getval(fname, 'IDCTAB'))
    except KeyError:
        return False
    if idctab == 'N/A' or idctab == "":
        return False
    if os.path.exists(idctab):
        return True
    else:
        return False

def applyTDDCorr(fname, utddcorr):
    """
    The default value of tddcorr for all ACS images is True.
    This correction will be performed if all conditions below are True:
    - the user did not turn it off on the command line
    - the detector is WFC
    - the idc table specified in the primary header is available.
    """

    phdr = pyfits.getheader(fname)
    instrument = phdr['INSTRUME']
    try:
        detector = phdr['DETECTOR']
    except KeyError:
        detector = None
    try:
        tddswitch = phdr['TDDCORR']
    except KeyError:
        tddswitch = 'PERFORM'

    if instrument == 'ACS' and detector == 'WFC' and utddcorr == True and tddswitch == 'PERFORM':
        tddcorr = True
        try:
            idctab = phdr['IDCTAB']
        except KeyError:
            tddcorr = False
            #print "***IDCTAB keyword not found - not applying TDD correction***\n"
        if os.path.exists(fileutil.osfn(idctab)):
            tddcorr = True
        else:
            tddcorr = False
            #print "***IDCTAB file not found - not applying TDD correction***\n"
    else:
        tddcorr = False

    return tddcorr

def applyNpolCorr(fname, unpolcorr):
    """
    Determines whether non-polynomial distortion lookup tables should be added
    as extensions to the science file based on the 'NPOLFILE' keyword in the
    primary header and NPOLEXT kw in the first extension.
    This is a default correction and will always run in the pipeline.
    The file used to generate the extensions is
    recorded in the NPOLEXT keyword in the first science extension.
    If 'NPOLFILE' in the primary header is different from 'NPOLEXT' in the
    extension header and the file exists on disk and is a 'new type' npolfile,
    then the lookup tables will be updated as 'WCSDVARR' extensions.
    """
    applyNPOLCorr = True
    try:
        # get NPOLFILE kw from primary header
        fnpol0 = pyfits.getval(fname, 'NPOLFILE')
        if fnpol0 == 'N/A':
            return False
        fnpol0 = fileutil.osfn(fnpol0)
        if not fileutil.findFile(fnpol0):
            msg =  """\n\tKw "NPOLFILE" exists in primary header but file %s not found
                      Non-polynomial distortion correction will not be applied\n
                    """ % fnpol0
            logger.critical(msg)
            applyNPOLCorr = False
            return applyNPOLCorr
        try:
            # get NPOLEXT kw from first extension header
            fnpol1 = pyfits.getval(fname, 'NPOLEXT', ext=1)
            fnpol1 = fileutil.osfn(fnpol1)
            if fnpol1 and fileutil.findFile(fnpol1):
                if fnpol0 != fnpol1:
                    applyNPOLCorr = True
                else:
                    msg = """\n\tNPOLEXT with the same value as NPOLFILE found in first extension.
                             NPOL correction will not be applied."""
                    logger.info(msg)
                    applyNPOLCorr = False
            else:
                # npl file defined in first extension may not be found
                # but if a valid kw exists in the primary header, non-polynomial
                #distortion correction should be applied.
                applyNPOLCorr = True
        except KeyError:
            # the case of "NPOLFILE" kw present in primary header but "NPOLEXT" missing
            # in first extension header
            applyNPOLCorr = True
    except KeyError:
        logger.info('\n\t"NPOLFILE" keyword not found in primary header')
        applyNPOLCorr = False
        return applyNPOLCorr

    if isOldStyleDGEO(fname, fnpol0):
            applyNPOLCorr = False
    return (applyNPOLCorr and unpolcorr)

def isOldStyleDGEO(fname, dgname):
    # checks if the file defined in a NPOLFILE kw is a full size
    # (old style) image

    sci_hdr = pyfits.getheader(fname, ext=1)
    dgeo_hdr = pyfits.getheader(dgname, ext=1)
    sci_naxis1 = sci_hdr['NAXIS1']
    sci_naxis2 = sci_hdr['NAXIS2']
    dg_naxis1 = dgeo_hdr['NAXIS1']
    dg_naxis2 = dgeo_hdr['NAXIS2']
    if sci_naxis1 <= dg_naxis1 or sci_naxis2 <= dg_naxis2:
        msg = """\n\tOnly full size (old style) DGEO file was found.\n
                 Non-polynomial distortion  correction will not be applied."""
        logger.critical(msg)
        return True
    else:
        return False

def applyD2ImCorr(fname, d2imcorr):
    applyD2IMCorr = True
    try:
        # get D2IMFILE kw from primary header
        fd2im0 = pyfits.getval(fname, 'D2IMFILE')
        if fd2im0 == 'N/A':
            return False
        fd2im0 = fileutil.osfn(fd2im0)
        if not fileutil.findFile(fd2im0):
            msg = """\n\tKw D2IMFILE exists in primary header but file %s not found\n
                     Detector to image correction will not be applied\n""" % fd2im0
            logger.critical(msg)
            print msg
            applyD2IMCorr = False
            return applyD2IMCorr
        try:
            # get D2IMEXT kw from first extension header
            fd2imext = pyfits.getval(fname, 'D2IMEXT', ext=1)
            fd2imext = fileutil.osfn(fd2imext)
            if fd2imext and fileutil.findFile(fd2imext):
                if fd2im0 != fd2imext:
                    applyD2IMCorr = True
                else:
                    applyD2IMCorr = False
            else:
                # D2IM file defined in first extension may not be found
                # but if a valid kw exists in the primary header,
                # detector to image correction should be applied.
                applyD2IMCorr = True
        except KeyError:
            # the case of D2IMFILE kw present in primary header but D2IMEXT missing
            # in first extension header
            applyD2IMCorr = True
    except KeyError:
        print 'D2IMFILE keyword not found in primary header'
        applyD2IMCorr = False
        return applyD2IMCorr

