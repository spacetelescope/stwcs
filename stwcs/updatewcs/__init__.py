from __future__ import absolute_import, division, print_function

import atexit
from astropy.io import fits
from .. import wcsutil
#from ..wcsutil.hwstwcs import HSTWCS

from .. import __version__

from astropy import wcs as pywcs
import astropy
from astropy import log
default_log_level = log.getEffectiveLevel()

from . import utils, corrections
from . import npol, det2im
from stsci.tools import parseinput, fileutil
from . import apply_corrections

import time
import logging
logger = logging.getLogger('stwcs.updatewcs')

atexit.register(logging.shutdown)

# Note: The order of corrections is important


def updatewcs(input, vacorr=True, tddcorr=True, npolcorr=True, d2imcorr=True,
              checkfiles=True, verbose=False):
    """

    Updates HST science files with the best available calibration information.
    This allows users to retrieve from the archive self contained science files
    which do not require additional reference files.

    Basic WCS keywords are updated in the process and new keywords (following WCS
    Paper IV and the SIP convention) as well as new extensions are added to the science files.


    Example
    -------
    >>>from stwcs import updatewcs
    >>>updatewcs.updatewcs(filename)

    Dependencies
    ------------
    `stsci.tools`
    `astropy.io.fits`
    `astropy.wcs`

    Parameters
    ----------
    input: a python list of file names or a string (wild card characters allowed)
             input files may be in fits, geis or waiver fits format
    vacorr: boolean
              If True, vecocity aberration correction will be applied
    tddcorr: boolean
             If True, time dependent distortion correction will be applied
    npolcorr: boolean
              If True, a Lookup table distortion will be applied
    d2imcorr: boolean
              If True, detector to image correction will be applied
    checkfiles: boolean
              If True, the format of the input files will be checked,
              geis and waiver fits files will be converted to MEF format.
              Default value is True for standalone mode.
    """
    if not verbose:
        logger.setLevel(100)
    else:
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        log_filename = 'stwcs.log'
        fh = logging.FileHandler(log_filename, mode='w')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logger.setLevel(verbose)
    args = "vacorr=%s, tddcorr=%s, npolcorr=%s, d2imcorr=%s, checkfiles=%s, \
    " % (str(vacorr), str(tddcorr), str(npolcorr),
         str(d2imcorr), str(checkfiles))
    logger.info('\n\tStarting UPDATEWCS: %s', time.asctime())

    files = parseinput.parseinput(input)[0]
    logger.info("\n\tInput files: %s, " % [i for i in files])
    logger.info("\n\tInput arguments: %s" % args)
    if checkfiles:
        files = checkFiles(files)
        if not files:
            print('No valid input, quitting ...\n')
            return

    for f in files:
        acorr = apply_corrections.setCorrections(f, vacorr=vacorr, tddcorr=tddcorr,
                                                 npolcorr=npolcorr, d2imcorr=d2imcorr)
        if 'MakeWCS' in acorr and newIDCTAB(f):
            logger.warning("\n\tNew IDCTAB file detected. All current WCSs will be deleted")
            cleanWCS(f)

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
    logger.info("Allowed corrections: {0}".format(allowed_corr))
    f = fits.open(fname, mode='update')
    f.readall()
    # Determine the reference chip and create the reference HSTWCS object
    nrefchip, nrefext = getNrefchip(f)
    log.setLevel("WARNING")
    wcsutil.restoreWCS(f, nrefext, wcskey='O')
    log.setLevel(default_log_level)
    rwcs = wcsutil.HSTWCS(fobj=f, ext=nrefext)
    rwcs.readModel(update=True, header=f[nrefext].header)

    if 'DET2IMCorr' in allowed_corr:
        kw2update = det2im.DET2IMCorr.updateWCS(f)
        for kw in kw2update:
            f[1].header[kw] = kw2update[kw]

    for i in range(len(f))[1:]:
        extn = f[i]

        if 'extname' in extn.header:
            extname = extn.header['extname'].lower()
            if extname == 'sci':
                log.setLevel('WARNING')
                wcsutil.restoreWCS(f, ext=i, wcskey='O')
                log.setLevel(default_log_level)
                sciextver = extn.header['extver']
                ref_wcs = rwcs.deepcopy()
                hdr = extn.header
                ext_wcs = wcsutil.HSTWCS(fobj=f, ext=i)
                # check if it exists first!!!
                # 'O ' can be safely archived again because it has been restored first.
                wcsutil.archiveWCS(f, ext=i, wcskey="O", wcsname="OPUS", reusekey=True)
                ext_wcs.readModel(update=True, header=hdr)
                for c in allowed_corr:
                    if c != 'NPOLCorr' and c != 'DET2IMCorr':
                        corr_klass = corrections.__getattribute__(c)
                        kw2update = corr_klass.updateWCS(ext_wcs, ref_wcs)
                        for kw in kw2update:
                            hdr[kw] = kw2update[kw]
                # give the primary WCS a WCSNAME value
                idcname = f[0].header.get('IDCTAB', " ")
                if idcname.strip() and 'idc.fits' in idcname:
                    wname = ''.join(['IDC_',
                                     utils.extract_rootname(idcname, suffix='_idc')])
                else: wname = " "
                hdr['WCSNAME'] = wname

            elif extname in ['err', 'dq', 'sdq', 'samp', 'time']:
                cextver = extn.header['extver']
                if cextver == sciextver:
                    hdr = f[('SCI', sciextver)].header
                    w = pywcs.WCS(hdr, f)
                    copyWCS(w, extn.header)

            else:
                continue

    if 'NPOLCorr' in allowed_corr:
        kw2update = npol.NPOLCorr.updateWCS(f)
        for kw in kw2update:
            f[1].header[kw] = kw2update[kw]
    # Finally record the version of the software which updated the WCS
    if 'HISTORY' in f[0].header:
        f[0].header.set('UPWCSVER', value=__version__,
                        comment="Version of STWCS used to updated the WCS",
                        before='HISTORY')
        f[0].header.set('PYWCSVER', value=astropy.__version__,
                        comment="Version of PYWCS used to updated the WCS",
                        before='HISTORY')
    elif 'ASN_MTYP' in f[0].header:
        f[0].header.set('UPWCSVER', value=__version__,
                        comment="Version of STWCS used to updated the WCS",
                        after='ASN_MTYP')
        f[0].header.set('PYWCSVER', value=astropy.__version__,
                        comment="Version of PYWCS used to updated the WCS",
                        after='ASN_MTYP')
    else:
        # Find index of last non-blank card, and insert this new keyword after that card
        for i in range(len(f[0].header) - 1, 0, -1):
            if f[0].header[i].strip() != '':
                break
            f[0].header.set('UPWCSVER', __version__,
                            "Version of STWCS used to updated the WCS",
                            after=i)
            f[0].header.set('PYWCSVER', astropy.__version__,
                            "Version of PYWCS used to updated the WCS",
                            after=i)
    # add additional keywords to be used by headerlets
    distdict = utils.construct_distname(f, rwcs)
    f[0].header['DISTNAME'] = distdict['DISTNAME']
    f[0].header['SIPNAME'] = distdict['SIPNAME']
    # Make sure NEXTEND keyword remains accurate
    f[0].header['NEXTEND'] = len(f) - 1
    f.close()


def copyWCS(w, ehdr):
    """
    This is a convenience function to copy a WCS object
    to a header as a primary WCS. It is used only to copy the
    WCS of the 'SCI' extension to the headers of 'ERR', 'DQ', 'SDQ',
    'TIME' or 'SAMP' extensions.
    """
    log.setLevel('WARNING')
    hwcs = w.to_header()
    log.setLevel(default_log_level)

    if w.wcs.has_cd():
        wcsutil.pc2cd(hwcs)
    for k in hwcs.keys():
        key = k[:7]
        ehdr[key] = hwcs[k]


def getNrefchip(fobj):
    """

    Finds which FITS extension holds the reference chip.

    The reference chip has EXTNAME='SCI', can be in any extension and
    is instrument specific. This functions provides mappings between
    sci extensions, chips and fits extensions.
    In the case of a subarray when the reference chip is missing, the
    first 'SCI' extension is the reference chip.

    Parameters
    ----------
    fobj: `astropy.io.fits.HDUList` object
    """
    nrefext = 1
    nrefchip = 1
    instrument = fobj[0].header['INSTRUME']

    if instrument == 'WFPC2':
        chipkw = 'DETECTOR'
        extvers = [("SCI", img.header['EXTVER']) for img in
                   fobj[1:] if img.header['EXTNAME'].lower() == 'sci']
        detectors = [img.header[chipkw] for img in fobj[1:] if
                     img.header['EXTNAME'].lower() == 'sci']
        fitsext = [i for i in range(len(fobj))[1:] if
                   fobj[i].header['EXTNAME'].lower() == 'sci']
        det2ext = dict(list(zip(detectors, extvers)))
        ext2det = dict(list(zip(extvers, detectors)))
        ext2fitsext = dict(list(zip(extvers, fitsext)))

        if 3 not in detectors:
            nrefchip = ext2det.pop(extvers[0])
            nrefext = ext2fitsext.pop(extvers[0])
        else:
            nrefchip = 3
            extname = det2ext.pop(nrefchip)
            nrefext = ext2fitsext.pop(extname)

    elif (instrument == 'ACS' and fobj[0].header['DETECTOR'] == 'WFC') or \
         (instrument == 'WFC3' and fobj[0].header['DETECTOR'] == 'UVIS'):
        chipkw = 'CCDCHIP'
        extvers = [("SCI", img.header['EXTVER']) for img in
                   fobj[1:] if img.header['EXTNAME'].lower() == 'sci']
        detectors = [img.header[chipkw] for img in fobj[1:] if
                     img.header['EXTNAME'].lower() == 'sci']
        fitsext = [i for i in range(len(fobj))[1:] if
                   fobj[i].header['EXTNAME'].lower() == 'sci']
        det2ext = dict(list(zip(detectors, extvers)))
        ext2det = dict(list(zip(extvers, detectors)))
        ext2fitsext = dict(list(zip(extvers, fitsext)))

        if 2 not in detectors:
            nrefchip = ext2det.pop(extvers[0])
            nrefext = ext2fitsext.pop(extvers[0])
        else:
            nrefchip = 2
            extname = det2ext.pop(nrefchip)
            nrefext = ext2fitsext.pop(extname)
    else:
        for i in range(len(fobj)):
            extname = fobj[i].header.get('EXTNAME', None)
            if extname is not None and extname.lower == 'sci':
                nrefext = i
                break

    return nrefchip, nrefext


def checkFiles(input):
    """
    Checks that input files are in the correct format.
    Converts geis and waiver fits files to multiextension fits.
    """
    from stsci.tools.check_files import geis2mef, waiver2mef, checkFiles
    logger.info("\n\tChecking files %s" % input)
    removed_files = []
    newfiles = []
    if not isinstance(input, list):
        input = [input]

    for file in input:
        try:
                imgfits, imgtype = fileutil.isFits(file)
        except IOError:
            logger.warning("File {0} could not be found, removing it from the"
                           "input list.".format(file))
            removed_files.append(file)
            continue
        # Check for existence of waiver FITS input, and quit if found.
        # Or should we print a warning and continue but not use that file
        if imgfits:
            if imgtype == 'waiver':
                newfilename = waiver2mef(file, convert_dq=True)
                if newfilename is None:
                    logger.warning("Could not convert waiver to mef."
                                   "Removing file {0} from input list".format(file))
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
            if newfilename is None:
                logger.warning("Could not convert file {0} from geis to mef, removing it from input list".format(file))
                removed_files.append(file)
            else:
                newfiles.append(newfilename)
    if removed_files:
        logger.warning('\n\tThe following files will be removed from the list of files',
                       'to be processed %s' % removed_files)

    newfiles = checkFiles(newfiles)[0]
    logger.info("\n\tThese files passed the input check and will be processed: %s" % newfiles)
    return newfiles


def newIDCTAB(fname):
    # When this is called we know there's a kw IDCTAB in the header
    hdul = fits.open(fname)
    idctab = fileutil.osfn(hdul[0].header['IDCTAB'])
    try:
        # check for the presence of IDCTAB in the first extension
        oldidctab = fileutil.osfn(hdul[1].header['IDCTAB'])
    except KeyError:
        return False
    if idctab == oldidctab:
        return False
    else:
        return True


def cleanWCS(fname):
    # A new IDCTAB means all previously computed WCS's are invalid
    # We are deleting all of them except the original OPUS WCS.
    f = fits.open(fname, mode='update')
    keys = wcsutil.wcskeys(f[1].header)
    # Remove the primary WCS from the list
    try:
        keys.remove(' ')
    except ValueError:
        pass
    fext = list(range(1, len(f)))
    for key in keys:
        try:
            wcsutil.deleteWCS(fname, ext=fext, wcskey=key)
        except KeyError:
            # Some extensions don't have the alternate (or any) WCS keywords
            continue


def getCorrections(instrument):
    """
    Print corrections available for an instrument

    :Parameters:
    `instrument`: string, one of 'WFPC2', 'NICMOS', 'STIS', 'ACS', 'WFC3'
    """
    acorr = apply_corrections.allowed_corrections[instrument]

    print("The following corrections will be performed for instrument %s\n" % instrument)
    for c in acorr: print(c, ': ', apply_corrections.cnames[c])
