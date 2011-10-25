from __future__ import division
import logging
import os
import string
import textwrap
import copy
import tarfile
import tempfile
import time
import warnings
from cStringIO import StringIO

import numpy as np
import pyfits

import altwcs
import wcscorr
from hstwcs import HSTWCS
from mappings import basic_wcs
from stwcs.updatewcs import utils

from stsci.tools.fileutil import countExtn
from stsci.tools import fileutil as fu

#### Logging support functions
module_logger = logging.getLogger('headerlet')

import atexit
atexit.register(logging.shutdown)

FITS_STD_KW = ['XTENSION', 'BITPIX', 'NAXIS', 'PCOUNT',
             'GCOUNT','EXTNAME', 'EXTVER', 'ORIGIN',
             'INHERIT', 'DATE', 'IRAF-TLM']

DEFAULT_SUMMARY_COLS = ['HDRNAME','WCSNAME','DISTNAME','AUTHOR','DATE',
                        'SIPNAME','NPOLFILE','D2IMFILE','DESCRIP']
COLUMN_DICT = {'vals':[],'width':[]}
COLUMN_FMT = '{:<{width}}'

def setLogger(logger, level, mode='w'):
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    log_filename = 'headerlet.log'
    fh = logging.FileHandler(log_filename, mode=mode)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.setLevel(level)
    
def initLogging(function_name, logger = None, level=100, verbose=False, logmode='w'):
    """ Initialize logging for a function 
    
    Parameters
    ----------
    function_name: string
        Name of function which will be recorded in log
    level: int
        Logging level
    verbose: bool
        If True, will set logging to report more activity
    """
    if logger is None:
        logger = module_logger
        
    if verbose:
        setLogger(logger, verbose, mode=logmode)
    else:
        logger.setLevel(level)

    logger.info("Starting %s: %s" % (function_name, time.asctime()))
    

#### Utility functions
def is_par_blank(par):
    return par in ['',' ','INDEF',"None",None]

def parseFilename(fname,mode='readonly'):
    """
    Interprets the input as either a filename of a file that needs to be opened
    or a PyFITS object.
    
    Parameters
    ----------
    fname: string, pyfits.HDUList
        Input pointing to a file or PyFITS object. An input filename (str) will
        be expanded as necessary to interpret any environmental variables 
        included in the filename.
    mode: string
        Specifies what PyFITS mode to use when opening the file, if it needs
        to open the file at all [Default: 'readonly']
        
    Returns
    -------
    fobj: pyfits.HDUList
        PyFITS handle for input file
    fname: string
        Name of input file
    close_fobj: bool
        Flag specifying whether or not fobj needs to be closed since it was 
        opened by this function. This allows a program to know whether they 
        need to worry about closing the PyFITS object as opposed to letting
        the higher level interface close the object.
    """
    close_fobj = False
    if not isinstance(fname,list):
        if isinstance(fname,str):
            fname = fu.osfn(fname)
        fobj = pyfits.open(fname,mode=mode)
        close_fobj = True
    else:
        fobj = fname
        if hasattr(fobj,'filename'):
            fname = fobj.filename()
        else:
            fname = ''
    return fobj,fname,close_fobj

def getHeaderletKwNames(fobj,kw='HDRNAME'):
    """
    Returns a list of specified keywords from all HeaderletHDU 
    extensions in a science file. 
    
    Parameters
    ----------
    fobj: string, pyfits.HDUList
    kw: str
        Name of keyword to be read and reported
    """
    fobj,fname,open_fobj = parseFilename(fobj)

    hdrnames = []
    for ext in fobj:
        if isinstance(ext,pyfits.hdu.base.NonstandardExtHDU):
            hdrnames.append(ext.header[kw])

    if open_fobj: 
        fobj.close()

    return hdrnames



def findHeaderletHDUs(fobj, hdrext=None, hdrname=None, distname=None, strict=True):
    """ 
    Returns all HeaderletHDU extensions in a science file that matches
    the inputs specified by the user.  If no hdrext, hdrname or distname are
    specified, this function will return a list of all HeaderletHDU objects.
    
    Parameters
    ----------
    fobj: string, pyfits.HDUList
        Name of FITS file or open pyfits object (pyfits.HDUList instance)
    hdrext: int, tuple or None
        index number(EXTVER) or extension tuple of HeaderletHDU to be returned
    hdrname: string
        value of HDRNAME for HeaderletHDU to be returned
    distname: string
        value of DISTNAME for HeaderletHDUs to be returned
    strict: bool [Default: True]
        Specifies whether or not at least one parameter needs to be provided
        If False, all extension indices returned if hdrext, hdrname and distname
        are all None. If True and hdrext, hdrname, and distname are all None,
        raise an Exception requiring one to be specified. 
    
    Returns
    -------
    hdrlets: list
        A list of all matching HeaderletHDU extension indices (could be just one)
    """    
    get_all = False
    if hdrext is None and hdrname is None and distname is None:
        if not strict:
            get_all = True
        else:
            print '====================================================='
            print 'No valid Headerlet extension specified.'
            print '  Either "hdrname", "hdrext", or "distname" needs to be specified.'
            print '====================================================='
            raise ValueError
            
    fobj,fname,open_fobj = parseFilename(fobj)

    hdrlets = []
    if hdrext is not None and isinstance(hdrext,int):
        if hdrext in range(len(fobj)): # insure specified hdrext is in fobj
            if isinstance(fobj[hdrext],pyfits.hdu.base.NonstandardExtHDU) and \
                fobj[hdrext].header['EXTNAME'] == 'HDRLET':
                    hdrlets.append(hdrext)
    else:
        for ext in fobj:            
            if isinstance(ext,pyfits.hdu.base.NonstandardExtHDU):
                if get_all:
                    hdrlets.append(fobj.index(ext))
                else:
                    if hdrext is not None:
                        if isinstance(hdrext,tuple):
                            hdrextname = hdrext[0]
                            hdrextnum = hdrext[1]
                        else:
                            hdrextname = 'HDRLET'
                            hdrextnum = hdrext
                    hdrext_match = ((hdrext is not None) and 
                                    (hdrextnum == ext.header['EXTVER']) and
                                    (hdrextname == ext.header['EXTNAME']))
                    hdrname_match = ((hdrname is not None) and 
                                        (hdrname == ext.header['HDRNAME']))
                    distname_match = ((distname is not None) and 
                                        (distname == ext.header['DISTNAME']))
                    if hdrext_match or hdrname_match or distname_match:
                        hdrlets.append(fobj.index(ext))
                
    if open_fobj: 
        fobj.close()
        
    if len(hdrlets) == 0:
        if hdrname:
            kwerr = 'hdrname'
            kwval = hdrname
        elif hdrext:
            kwerr = 'hdrext'
            kwval = hdrext
        else:
            kwerr = 'distname'
            kwval = distname
        print '====================================================='
        print 'No valid Headerlet extension found!'
        print '  "%s" = %s not found in %s.'%(kwerr,kwval,fname)
        print '====================================================='
        raise ValueError

    return hdrlets

def verifyHdrnameIsUnique(fobj,hdrname):
    """ 
    Verifies that no other HeaderletHDU extension has the specified hdrname.
    
    Parameters
    ----------
    fobj: string, pyfits.HDUList
        Name of FITS file or open pyfits object (pyfits.HDUList instance)
    hdrname: string
        value of HDRNAME for HeaderletHDU to be compared as unique
        
    Returns
    -------
    unique: bool
        If True, no other HeaderletHDU has the specified HDRNAME value
    """
    hdrnames_list = getHeaderletKwNames(fobj)
    unique = not(hdrname in hdrnames_list)
    
    return unique

    
def isWCSIdentical(scifile, file2, verbose=False):
    """
    Compares the WCS solution of 2 files.

    Parameters
    ----------
    scifile: file1
    file2: file2
    verbose: False or a python logging level
             (one of 'INFO', 'DEBUG' logging levels)
             (an integer representing a logging level)

    Notes
    -----
    These can be 2 science observations or 2 headerlets
    or a science observation and a headerlet. The two files
    have the same WCS solution if the following are the same:

    - rootname/destim
    - primary WCS
    - SIP coefficients
    - NPOL distortion
    - D2IM correction
    - Velocity aberation

    """
    initLogging('isWCSIdentical', level=100, verbose=False)

    result = True
    numsci1 = max(countExtn(scifile), countExtn(scifile, 'SIPWCS'))
    numsci2 = max(countExtn(file2), countExtn(file2, 'SIPWCS'))

    if numsci1 == 0 or numsci2 == 0 or numsci1 != numsci2:
        module_logger.info("Number of SCI and SIPWCS extensions do not match.")
        result = False

    if getRootname(scifile) != getRootname(file2):
        module_logger.info('Rootnames do not match.')
        result = False
    try:
        extname1 = pyfits.getval(scifile, 'EXTNAME', ext=('SCI', 1))
    except KeyError:
        extname1 = 'SIPWCS'
    try:
        extname2 = pyfits.getval(file2, 'EXTNAME', ext=('SCI', 1))
    except KeyError:
        extname2 = 'SIPWCS'

    for i in range(1, numsci1 + 1):
        w1 = HSTWCS(scifile, ext=(extname1, i))
        w2 = HSTWCS(file2, ext=(extname2, i))
        if not np.allclose(w1.wcs.crval, w2.wcs.crval, rtol=1e-7) or \
        not np.allclose(w1.wcs.crpix, w2.wcs.crpix, rtol=1e-7)  or \
        not np.allclose(w1.wcs.cd, w2.wcs.cd, rtol=1e-7) or \
        not (np.array(w1.wcs.ctype) == np.array(w2.wcs.ctype)).all():
            module_logger.info('Primary WCSs do not match')
            result = False
        if w1.sip or w2.sip:
            if (w2.sip and not w1.sip) or (w1.sip and not w2.sip) or \
               not np.allclose(w1.sip.a, w2.sip.a, rtol=1e-7) or \
               not np.allclose(w1.sip.b, w2.sip.b, rtol=1e-7):
                module_logger.info('SIP coefficients do not match')
                result = False
        if w1.cpdis1 or w2.cpdis1:
            if w1.cpdis1 and not w2.cpdis1 or \
                w2.cpdis1 and not w1.cpdis1 or \
                not np.allclose(w1.cpdis1.data, w2.cpdis1.data):
                module_logger.info('NPOL distortions do not match')
                result = False
        if w1.cpdis2 or w2.cpdis2:
            if w1.cpdis2 and not w2.cpdis2 or \
                w2.cpdis2 and not w1.cpdis2 or \
                not np.allclose(w1.cpdis2.data, w2.cpdis2.data):
                module_logger.info('NPOL distortions do not match')
                result = False
        if w1.det2im1 or w2.det2im1:
            if w1.det2im1 and not w2.det2im1 or \
                w2.det2im1 and not w1.det2im1 or\
                not np.allclose(w1.det2im1.data, w2.det2im1.data):
                module_logger.info('Det2Im corrections do not match')
                result =  False
        if w1.det2im2 or w2.det2im2:
            if w1.det2im2 and not w2.det2im2 or \
                w2.det2im2 and not w1.det2im2 or\
                not np.allclose(w1.det2im2.data, w2.det2im2.data):
                module_logger.info('Det2Im corrections do not match')
                result = False
        if w1.vafactor != w2.vafactor:
            module_logger.info('VA factors do not match')
            result = False

    return result

def updateRefFiles(source, dest, verbose=False):
    """
    Update the reference files name in the primary header of 'dest'
    using values from 'source'

    Parameters
    ----------
    source: pyfits.Header.ascardlist
    dest:   pyfits.Header.ascardlist
    """
    module_logger.info("Updating reference files")
    phdukw = {'IDCTAB': True,
            'NPOLFILE': True,
            'D2IMFILE': True}

    try:
        wind = dest.index_of('HISTORY')
    except KeyError:
        wind = len(dest)
    for key in phdukw.keys():
        try:
            value = source[key]
            dest.insert(wind, value)
        except KeyError:
            # TODO: I don't understand what the point of this is.  Is it meant
            # for logging purposes?  Right now it isn't used.
            phdukw[key] = False
    return phdukw

def getRootname(fname):
    """
    returns the value of ROOTNAME or DESTIM
    """

    try:
        rootname = pyfits.getval(fname, 'ROOTNAME')
    except KeyError:
        rootname = pyfits.getval(fname, 'DESTIM')
    return rootname

def mapFitsExt2HDUListInd(fname, extname):
    """
    Map FITS extensions with 'EXTNAME' to HDUList indexes.
    """

    f,fname,close_fobj = parseFilename(fname)

    d = {}
    for hdu in f:
        if 'EXTNAME' in hdu.header and hdu.header['EXTNAME'] == extname:
            extver = hdu.header['EXTVER']
            d[(extname, extver)] = f.index_of((extname, extver))
    if close_fobj:
        f.close()
    return d

def getHeaderletObj(fname, sciext='SCI', wcskey=' ',wcsname=None,hdrname=None, 
                    sipname=None, npolfile=None, d2imfile=None, 
                    author=None, descrip=None, history=None,
                    rms_ra = None, rms_dec = None, nmatch=None, catalog=None,
                    hdrletnum=1, verbose=100):
    """
    Generate and return a HeaderletHDU with EXTVER and HDRNAME set 
    
    """
    fobj,fname,open_fobj = parseFilename(fname)

    # Translate 'wcskey' value for PRIMARY WCS to valid altwcs value of ' '
    if wcskey == 'PRIMARY': wcskey = ' '
        
    hlt = create_headerlet(fobj, sciext=sciext, 
                             wcskey=wcskey, wcsname=wcsname, hdrname=hdrname, 
                             sipname=sipname, npolfile=npolfile, d2imfile=d2imfile, 
                             author=author, descrip=descrip, history=history,
                             rms_ra=rms_ra, rms_dec=rms_dec, nmatch=nmatch, 
                             catalog=catalog, verbose=verbose, logmode='a')    
    if open_fobj:
        fobj.close()

    return hlt

def print_summary(summary_cols, summary_dict, pad=2, maxwidth=None, idcol=None, 
                    output=None, clobber=True, quiet=False ):
    """ 
    Print out summary dictionary to STDOUT, and possibly an output file
    
    """
    nrows = None
    if idcol:
        nrows = len(idcol['vals'])

    # Find max width of each column
    column_widths = {}
    for kw in summary_dict:
        colwidth = np.array(summary_dict[kw]['width']).max()
        if maxwidth:
            colwidth = min(colwidth,maxwidth)
        column_widths[kw] = colwidth + pad
        if nrows is None:
            nrows = len(summary_dict[kw]['vals'])

    # print rows now
    outstr = ''
    # Start with column names
    if idcol:
        outstr += COLUMN_FMT.format(idcol['name'],width=idcol['width']+pad)
    for kw in summary_cols:
        outstr += COLUMN_FMT.format(kw,width=column_widths[kw])
    outstr += '\n'
    # Now, add a row for each headerlet
    for row in range(nrows):
        if idcol:
            outstr += COLUMN_FMT.format(idcol['vals'][row],width=idcol['width']+pad)
        for kw in summary_cols:
            val = summary_dict[kw]['vals'][row][:(column_widths[kw]-pad)]
            outstr += COLUMN_FMT.format(val, width=column_widths[kw])
        outstr += '\n'
    if not quiet:
        print outstr

    # If specified, write info to separate text file
    write_file = False
    if output:
        output = fu.osfn(output) # Expand any environment variables in filename
        write_file = True
        if os.path.exists(output):
            if clobber: 
                os.remove(output)
            else: 
                print 'WARNING: Not writing results to file!'
                print '         Output text file ',output,' already exists.'
                print '         Set "clobber" to True or move file before trying again.'
                write_file = False
        if write_file:
            fout = open(output,mode='w')
            fout.write(outstr)
            fout.close()    

#### Private utility functions
def _createPrimaryHDU(destim, hdrname, distname, wcsname, 
                            sipname, npolfile, d2imfile, upwcsver, pywcsver,
                            author, descrip, history):
    if author is None: author = ''
    if descrip is None: descrip = ''
    if history is None: history = ''
    phdu = pyfits.PrimaryHDU()
    phdu.header.update('DESTIM', destim,
                       comment='Destination observation root name')
    phdu.header.update('HDRNAME', hdrname, comment='Headerlet name')
    fmt="%Y-%m-%dT%H:%M:%S"
    phdu.header.update('DATE', time.strftime(fmt),
                       comment='Date FITS file was generated')
    phdu.header.update('WCSNAME', wcsname, comment='WCS name')
    phdu.header.update('DISTNAME', distname, comment='Distortion model name')
    phdu.header.update('SIPNAME', sipname, comment='origin of SIP polynomial distortion model')
    phdu.header.update('NPOLFILE', npolfile, comment='origin of non-polynmial distortion model')
    phdu.header.update('D2IMFILE', d2imfile, comment='origin of detector to image correction')
    phdu.header.update('AUTHOR', author, comment='headerlet created by this user')
    phdu.header.update('DESCRIP', descrip, comment='Short description of headerlet solution')
    
    # clean up history string in order to remove whitespace characters that
    # would cause problems with FITS
    if isinstance(history,list):
        history_str = ''
        for line in history:
            history_str += line
    else:
        history_str = history
    history_lines = textwrap.wrap(history_str,width=70)
    for hline in history_lines:
        phdu.header.add_history(hline)
    
    phdu.header.ascard.append(upwcsver)
    phdu.header.ascard.append(pywcsver)
    return phdu

    
#### Public Interface functions
def extract_headerlet(filename, output, extnum=None, hdrname=None,  
                    clobber=False, verbose=100):
    """
    Finds a headerlet extension in a science file and writes it out as 
    a headerlet FITS file.
     
    If both hdrname and extnum are given they should match, if not 
    raise an Exception
    
    Parameters
    ----------
    filename: string or HDUList
           Either a filename or PyFITS HDUList object for the input science file
            An input filename (str) will be expanded as necessary to interpret 
            any environmental variables included in the filename.
    output: string
           Filename or just rootname of output headerlet FITS file
           If string does not contain '.fits', it will create a filename with
           '_hlet.fits' suffix
    extnum: int
           Extension number which contains the headerlet to be written out
    hdrname: string
           Unique name for headerlet, stored as the HDRNAME keyword
           It stops if a value is not provided and no extnum has been specified  
    clobber: bool
        If output file already exists, this parameter specifies whether or not
        to overwrite that file [Default: False]
    verbose: int
             python logging level
            
    """
    
    initLogging('extract_headerlet', verbose=verbose)
    fobj,fname,close_fobj = parseFilename(filename)

    if hdrname in ['',' ',None, 'INDEF']:
        if  extnum is None:
            print 'No valid headerlet specified! Quitting...'
            if close_fobj: fobj.close()
        else:
            hdrhdu = fobj[extnum]
    else:
        extnum = findHeaderletHDUs(fobj,hdrname=hdrname)[0]
        hdrhdu = fobj[extnum]
    hdrlet = hdrhdu.headerlet
    
    if '.fits' in output:
        outname = output
    else:
        outname = '%s_hlet.fits'%(output)
    hdrlet.write_to(outname)

    if close_fobj:
        fobj.close()

def write_headerlet(filename, hdrname, output=None, sciext='SCI', 
                        wcsname=None, wcskey=None, destim=None,
                        sipname=None, npolfile=None, d2imfile=None, 
                        author=None, descrip=None, history=None,
                        rms_ra=None, rms_dec=None, nmatch=None, catalog=None,
                        attach=True, clobber=False):

    """
    Save a WCS as a headerlet FITS file.

    This function will create a headerlet, write out the headerlet to a 
    separate headerlet file, then, optionally, attach it as an extension 
    to the science image (if it has not already been archived)     
    
    Either wcsname or wcskey must be provided; if both are given, they must 
    match a valid WCS.

    Updates wcscorr if necessary.

    Parameters
    ----------
    filename: string or HDUList
           Either a filename or PyFITS HDUList object for the input science file
            An input filename (str) will be expanded as necessary to interpret 
            any environmental variables included in the filename.
    hdrname: string
        Unique name for this headerlet, stored as HDRNAME keyword 
    output: string or None
        Filename or just rootname of output headerlet FITS file
        If string does not contain '.fits', it will create a filename with
        '_hlet.fits' suffix
        If None, a default filename based on the input filename will be 
        generated for the headerlet FITS filename
    sciext: string
        name (EXTNAME) of extension that contains WCS to be saved
    wcsname: string
        name of WCS to be archived, if " ": stop
    wcskey: one of A...Z or " " or "PRIMARY" 
        if " " or "PRIMARY" - archive the primary WCS
    destim: string
        DESTIM keyword
        if  NOne, use ROOTNAME or science file name
    sipname: string or None (default)
         Name of unique file where the polynomial distortion coefficients were
         read from. If None, the behavior is:
         The code looks for a keyword 'SIPNAME' in the science header
         If not found, for HST it defaults to 'IDCTAB'
         If there is no SIP model the value is 'NOMODEL'
         If there is a SIP model but no SIPNAME, it is set to 'UNKNOWN'
    npolfile: string or None (default)
         Name of a unique file where the non-polynomial distortion was stored.
         If None:
         The code looks for 'NPOLFILE' in science header.
         If 'NPOLFILE' was not found and there is no npol model, it is set to 'NOMODEL'
         If npol model exists, it is set to 'UNKNOWN'
    d2imfile: string
         Name of a unique file where the detector to image correction was
         stored. If None:
         The code looks for 'D2IMFILE' in the science header.
         If 'D2IMFILE' is not found and there is no d2im correction,
         it is set to 'NOMODEL'
         If d2im correction exists, but 'D2IMFILE' is missing from science
         header, it is set to 'UNKNOWN'
    author: string
        Name of user who created the headerlet, added as 'AUTHOR' keyword
        to headerlet PRIMARY header
    descrip: string
        Short description of the solution provided by the headerlet
        This description will be added as the single 'DESCRIP' keyword
        to the headerlet PRIMARY header 
    history: filename, string or list of strings
        Long (possibly multi-line) description of the solution provided
        by the headerlet. These comments will be added as 'HISTORY' cards
        to the headerlet PRIMARY header
        If filename is specified, it will format and attach all text from
        that file as the history.
    attach: bool
        Specify whether or not to attach this headerlet as a new extension 
        It will verify that no other headerlet extension has been created with
        the same 'hdrname' value.
    clobber: bool
        If output file already exists, this parameter specifies whether or not
        to overwrite that file [Default: False]
    """
    initLogging('write_headerlet')

    if wcsname in [None,' ','','INDEF'] and wcskey is None:
        print '='*60
        print '[write_headerlet]'
        print 'No valid WCS found found in %s.'%fname
        print '   A valid value for either "wcsname" or "wcskey" '
        print '   needs to be specified. '
        print '='*60
        raise ValueError
    
    if hdrname in [None, ' ','']:
        print '='*60
        print '[write_headerlet]'
        print 'No valid name for this headerlet was provided for %s.'%fname
        print '   A valid value for "hdrname" needs to be specified. '
        print '='*60
        raise ValueError
        
    # Translate 'wcskey' value for PRIMARY WCS to valid altwcs value of ' '
    if wcskey == 'PRIMARY': wcskey = ' '
    
    if attach: umode = 'update'
    else: umode='readonly'
    
    fobj,fname,close_fobj = parseFilename(filename,mode=umode)

    # Insure that WCSCORR table has been created with all original
    # WCS's recorded prior to adding the headerlet WCS
    wcscorr.init_wcscorr(fobj)

    numhlt = countExtn(fobj, 'HDRLET')

    if wcsname is None:
        scihdr = fobj[sciext,1].header
        wcsname = scihdr['wcsname'+wcskey]

    hdrletobj = getHeaderletObj(fobj,sciext=sciext,
                                wcsname=wcsname, wcskey=wcskey,
                                hdrname=hdrname,
                                sipname=sipname, npolfile=npolfile, 
                                d2imfile=d2imfile, author=author, 
                                descrip=descrip, history=history,
                                rms_ra=rms_ra, rms_dec=rms_dec, nmatch=nmatch, 
                                catalog=catalog,
                                hdrletnum=numhlt + 1, verbose=False) 
            
    if attach:
        # Check to see whether or not a HeaderletHDU with this hdrname already 
        # exists
        hdrnames = getHeaderletKwNames(fobj)
        if hdrname not in hdrnames:
            hdrlet_hdu = HeaderletHDU.fromheaderlet(hdrletobj)

            if destim is not None:
                hdrlet_hdu[0].header['destim'] = destim

            fobj.append(hdrlet_hdu)
            
            # Update the WCSCORR table with new rows from the headerlet's WCSs
            wcscorr.update_wcscorr(fobj, hdrletobj, 'SIPWCS')

            fobj.flush()
        else:
            print 'WARNING:'
            print '    Headerlet with hdrname ',hdrname,' already archived for WCS ',wcsname
            print '    No new headerlet appended to ',fname,'.'
        
        
    if close_fobj:
        fobj.close()

    if output is None:
        # Generate default filename for headerlet FITS file
        output = fname[:fname.find('.fits')]
    if '.fits' not in output:
        output = output+'_hlet.fits'

    # If user specifies an output filename for headerlet, write it out
    if os.path.exists(output):
        if clobber:
            os.remove(output)
        else:
            print 'WARNING:'
            print '    Headerlet file ',output,' already written out!'
            print '    This headerlet file will not be created now.'

    hdrletobj.writeto(output)
    print 'Create Headerlet file: ',output
    
def create_headerlet(filename, sciext=None, hdrname=None, destim=None, wcskey=" ", wcsname=None, 
                     sipname=None, npolfile=None, d2imfile=None, 
                     author=None, descrip=None, history=None,
                     rms_ra=None, rms_dec = None, nmatch=None, catalog=None, 
                     verbose=100, logmode='w'):
    """
    Create a headerlet from a WCS in a science file   
    If both wcskey and wcsname are given they should match, if not 
    raise an Exception
    
    Parameters
    ----------
    filename: string or HDUList
           Either a filename or PyFITS HDUList object for the input science file
            An input filename (str) will be expanded as necessary to interpret 
            any environmental variables included in the filename.
    sciext: string or python list
           Extension in which the science data is. The headerlet will be created 
           from these extensions.
           If string - a valid EXTNAME is expected
           If list - a list of FITS extension numbers or extension tuples ('SCI', 1)
           is expected.
           If None, loops over all extensions in the file, including the primary [BROKEN-22Sept2011]
    hdrname: string
           value of HDRNAME keyword
           Takes the value from the HDRNAME<wcskey> keyword, if not available from WCSNAME<wcskey>
           It stops if neither is found in the science file and a value is not provided 
    destim: string or None
            name of file this headerlet can be applied to
            if None, use ROOTNAME keyword
    wcskey: char (A...Z) or " " or "PRIMARY" or None
            a char representing an alternate WCS to be used for the headerlet
            if " ", use the primary (default) 
            if None use wcsname           
    wcsname: string or None
            if wcskey is None use wcsname specified here to choose an alternate WCS for the headerlet
    sipname: string or None (default)
             Name of unique file where the polynomial distortion coefficients were
             read from. If None, the behavior is:
             The code looks for a keyword 'SIPNAME' in the science header
             If not found, for HST it defaults to 'IDCTAB'
             If there is no SIP model the value is 'NOMODEL'
             If there is a SIP model but no SIPNAME, it is set to 'UNKNOWN'
    npolfile: string or None (default)
             Name of a unique file where the non-polynomial distortion was stored.
             If None:
             The code looks for 'NPOLFILE' in science header.
             If 'NPOLFILE' was not found and there is no npol model, it is set to 'NOMODEL'
             If npol model exists, it is set to 'UNKNOWN'
    d2imfile: string
             Name of a unique file where the detector to image correction was
             stored. If None:
             The code looks for 'D2IMFILE' in the science header.
             If 'D2IMFILE' is not found and there is no d2im correction,
             it is set to 'NOMODEL'
             If d2im correction exists, but 'D2IMFILE' is missing from science
             header, it is set to 'UNKNOWN'
    author: string
            Name of user who created the headerlet, added as 'AUTHOR' keyword
            to headerlet PRIMARY header
    descrip: string
            Short description of the solution provided by the headerlet
            This description will be added as the single 'DESCRIP' keyword
            to the headerlet PRIMARY header 
    history: filename, string or list of strings
            Long (possibly multi-line) description of the solution provided
            by the headerlet. These comments will be added as 'HISTORY' cards
            to the headerlet PRIMARY header
            If filename is specified, it will format and attach all text from
            that file as the history.
    verbose: int
             python logging level
    logmode: 'w' or 'a'
             log file open mode
            
    Returns
    -------
    Headerlet object
    """
    
    initLogging('createHeaderlet', verbose=verbose)

    phdukw = {'IDCTAB': True,
            'NPOLFILE': True,
            'D2IMFILE': True}
    fobj,fname,close_file = parseFilename(filename)

    # Translate 'wcskey' value for PRIMARY WCS to valid altwcs value of ' '
    if wcskey == 'PRIMARY': wcskey = ' '
    
    # get all required keywords
    if destim is None:
        try:
            destim = fobj[0].header['ROOTNAME']
        except KeyError:
            destim = fname
            module_logger.info('DESTIM not provided')
            module_logger.info('Keyword "ROOTNAME" not found')
            module_logger.info('Using file name as DESTIM')
            
    if not hdrname:
        # check if HDRNAME<wcskey> is in header
        hdrname = "".join(["HDRNAME",wcskey.upper()])
        try:
            hdrname = fobj[1].header['HDRNAME']
        except KeyError:
            try:
                hdrname = fobj[1].header['WCSNAME'+wcskey]
                print 'Using default value for HDRNAME of "%d"',
                print ' derived from WCSNAME%d.'%(hdrname,wcskey)
            except KeyError, detail:
                message = "Required keywords 'HDRNAME' or 'WCSNAME' not found"
                module_logger.critical(message)
                print message, detail
    
    if not wcsname:
        wname = "".join(["WCSNAME",wcskey.upper()])
        if wname in fobj[1].header:
            wcsname = fobj[1].header[wname]
        else:
            message = "Missing required keyword 'WCSNAME'."
            module_logger.critical(message)
            print message, detail
            
    if not sipname:
        sipname = utils.build_sipname(fobj)
    
    if not npolfile:
        npolfile = utils.build_npolname(fobj)
            
    if not d2imfile:
        d2imfile = utils.build_d2imname(fobj)
    
    distname = utils.build_distname(sipname, npolfile, d2imfile)
    
    # get the version of STWCS used to create the WCS of the science file.
    try:
        upwcsver = fobj[0].header.ascard['STWCSVER']
    except KeyError:
        upwcsver = pyfits.Card("STWCSVER", " ",
                               "Version of STWCS used to update the WCS")
    try:
        pywcsver = fobj[0].header.ascard['PYWCSVER']
    except KeyError:
        pywcsver = pyfits.Card("PYWCSVER", " ",
                               "Version of PYWCS used to update the WCS")

    if not sciext:
        sciext = range(len(fobj))
    elif isinstance(sciext, str):
        numsciext = countExtn(fobj, sciext)
        sciext = [(sciext, i) for i in range(1, numsciext+1)]
    elif isinstance(sciext, list):
        pass
    else:
        raise ValueError("Expected sciext to be a list of FITS extensions with science data or a string of valid EXTNAME")

    if wcskey is 'O':
        message = "Warning: 'O' is a reserved key for the original WCS. Quitting..."
        module_logger.info(message)
        print message
        return

    # open file and parse comments
    if history not in ['',' ',None,'INDEF'] and os.path.isfile(history):
        f = open(fu.osfn(history))
        history = f.readlines()
        f.close()

    module_logger.debug("Data extensions form which to create headerlet:\n\t %s"
                 % (str(sciext)))
    hdul = pyfits.HDUList()
    phdu = _createPrimaryHDU(destim, hdrname, distname, wcsname, 
                             sipname, npolfile, d2imfile, upwcsver, pywcsver,
                             author, descrip, history)
    hdul.append(phdu)
    orient_comment = "positions angle of image y axis (deg. e of n)"
    
    if fu.isFits(fobj)[1] is not 'simple':
        
        for e in sciext:
            wkeys = altwcs.wcskeys(fname,ext=e)
            if wcskey != ' ' and wkeys > 0:
                if wcskey not in wkeys:
                    if verbose > 100:
                        module_logger.debug('No WCS with wcskey=%s found in extension %s.  Skipping...'%(wcskey,str(e)))
                    continue # skip any extension which does not have this wcskey
    
            # This reads in full model: alternate WCS keywords plus SIP
            hwcs = HSTWCS(fname,ext=e,wcskey=' ') 
            h = hwcs.wcs2header(sip2hdr=True)
            h.update('ORIENTAT',hwcs.orientat, comment=orient_comment)
            if wcskey != ' ':
                # Now read in specified linear WCS terms from alternate WCS
                try:
                    althdr = altwcs.convertAltWCS(fname,e,oldkey=wcskey,newkey=" ")
                    althdrwcs = HSTWCS(fname,e,wcskey=wcskey)
                    alt_orient = althdrwcs.orientat
                except KeyError:
                    continue # Skip over any extension which does not have a WCS
                althdr = althdr.ascard
                # Update full WCS with values from alternate WCS 
                for card in althdr:
                    h.update(card.key,card.value)
                h.update('ORIENTAT',alt_orientat, comment=orient_comment)
            h = h.ascard
            h.append(pyfits.Card(key='VAFACTOR', value=hwcs.vafactor,
                                 comment='Velocity aberration plate scale factor'))
            h.insert(0, pyfits.Card(key='EXTNAME', value='SIPWCS',
                                    comment='Extension name'))
            if isinstance(e, int): 
                if 'extver' in fobj[e].header:
                    val = fobj[e].header['extver']
                else:
                    val = e
            else: val = e[1]
            h.insert(1, pyfits.Card(key='EXTVER', value=val,
                                    comment='Extension version'))
            h.append(pyfits.Card("SCIEXT", str(e), 
                                 "Target science data extension"))
            fhdr = fobj[e].header.ascard
            if npolfile is not 'NOMODEL':
                cpdis = fhdr['CPDIS*...']
                for c in range(1, len(cpdis) + 1):
                    h.append(cpdis[c - 1])
                    dp = fhdr['DP%s*' % c]
                    h.extend(dp)
    
                    try:
                        h.append(fhdr['CPERROR%s' % c])
                    except KeyError:
                        pass
    
                try:
                    h.append(fhdr['NPOLEXT'])
                except KeyError:
                    pass
    
            if d2imfile is not 'NOMODEL':
                try:
                    h.append(fhdr['D2IMEXT'])
                except KeyError:
                    pass
    
                try:
                    h.append(fhdr['AXISCORR'])
                except KeyError:
                    print ("'D2IMFILE' kw exists but keyword 'AXISCORR' was not found in "
                                     "%s['SCI',%d]" % (fname, val))
                    module_logger.exception("'D2IMFILE' kw exists but keyword 'AXISCORR' was not found in "
                                     "%s['SCI',%d]" % (fname, val))
                    raise
    
                try:
                    h.append(fhdr['D2IMERR'])
                except KeyError:
                    h.append(pyfits.Card(key='DPERROR', value=0,
                                         comment='Maximum error of D2IMARR'))
    
            hdu = pyfits.ImageHDU(header=pyfits.Header(h))
            hdul.append(hdu)
    numwdvarr = countExtn(fname, 'WCSDVARR')
    numd2im = countExtn(fname, 'D2IMARR')
    for w in range(1, numwdvarr + 1):
        hdu = fobj[('WCSDVARR', w)].copy()
        hdul.append(hdu)
    for d in range(1, numd2im + 1):
        hdu = fobj[('D2IMARR', d)].copy()
        hdul.append(hdu)

    if close_file:
        fobj.close()
    return Headerlet(hdul,verbose=verbose, logmode='a')


def apply_headerlet_as_primary(filename, hdrlet, attach=True,archive=True,
                                force=False, verbose=False):
    """
    Apply headerlet 'hdrfile' to a science observation 'destfile' as the primary WCS

    Parameters
    ----------
    filename: string
             File name of science observation whose WCS solution will be updated
    hdrlet: string
             Headerlet file
    attach: boolean
            True (default): append headerlet to FITS file as a new extension.
    archive: boolean
            True (default): before updating, create a headerlet with the
            WCS old solution.
    force: boolean
            If True, this will cause the headerlet to replace the current PRIMARY
            WCS even if it has a different distortion model. [Default: False]
    verbose: False or a python logging level
             (one of 'INFO', 'DEBUG' logging levels)
             (an integer representing a logging level)
    """
    initLogging('apply_headerlet_as_primary', verbose=verbose)
    
    hlet = Headerlet(hdrlet, verbose=verbose, logmode='a')
    hlet.apply_as_primary(filename, attach=attach, archive=archive, force=force)

def apply_headerlet_as_alternate(filename, hdrlet, attach=True, 
                                wcskey=None, wcsname=None, verbose=False):
    """
    Apply headerlet to a science observation as an alternate WCS

    Parameters
    ----------
    filename: string
             File name of science observation whose WCS solution will be updated
    hdrlet: string
             Headerlet file
    attach: boolean
          flag indicating if the headerlet should be attached as a 
          HeaderletHDU to fobj. If True checks that HDRNAME is unique 
          in the fobj and stops if not.
    wcskey: string
          Key value (A-Z, except O) for this alternate WCS
          If None, the next available key will be used 
    wcsname: string
          Name to be assigned to this alternate WCS
          WCSNAME is a required keyword in a Headerlet but this allows the
          user to change it as desired.
    verbose: False or a python logging level
             (one of 'INFO', 'DEBUG' logging levels)
             (an integer representing a logging level)
    """
    initLogging('apply_headerlet_as_alternate', verbose=verbose)
    
    hlet = Headerlet(hdrlet, verbose=verbose, logmode='a')
    hlet.apply_as_alternate(filename, attach=attach, wcsname=wcsname, wcskey=wcskey)

def attach_headerlet(filename, hdrlet, verbose=False):
    """    
    Attach Headerlet as an HeaderletHDU to a science file   

    Parameters
    ----------
    filename: string, HDUList
            science file to which the headerlet should be applied
    hdrlet: string or Headerlet object
            string representing a headerlet file
    """
    initLogging('attach_headerlet', verbose=verbose)
    
    hlet = Headerlet(hdrlet, verbose=verbose, logmode='a')
    hlet.attach_to_file(filename)

def delete_headerlet(filename, hdrname=None, hdrext=None, distname=None):
    """
    Deletes HeaderletHDU(s) from a science file

    Notes
    -----
    One of hdrname, hdrext or distname should be given.
    If hdrname is given - delete a HeaderletHDU with a name HDRNAME from fobj.
    If hdrext is given - delete HeaderletHDU in extension.
    If distname is given - deletes all HeaderletHDUs with a specific distortion model from fobj.  
    Updates wcscorr

    Parameters
    ----------
    filename: string or HDUList
           Either a filename or PyFITS HDUList object for the input science file
            An input filename (str) will be expanded as necessary to interpret 
            any environmental variables included in the filename.
    hdrname: string or None
        HeaderletHDU primary header keyword HDRNAME
    hdrext: int, tuple or None
        HeaderletHDU FITS extension number
        tuple has the form ('HDRLET', 1)
    distname: string or None
        distortion model as specified in the DISTNAME keyword
    """    
    initLogging('delete_headerlet')
    hdrlet_ind = findHeaderletHDUs(filename,hdrname=hdrname, hdrext=hdrext, 
                            distname=distname)
    if len(hdrlet_ind) == 0:
        print 'ERROR: '
        print 'No HDUs deleted... No Headerlet HDUs found with '
        print '    hdrname = ',hdrname
        print '    hdrext  = ',hdrext
        print '   distname = ',distname
        print 'Please review input parameters and try again. '
        return

    fobj,fname,close_fobj = parseFilename(filename,mode='update')
    
    # delete row(s) from WCSCORR table now...
    #
    #
    if hdrname not in ['',' ',None,'INDEF']:
        selections = {'hdrname':hdrname}
    elif hdrname in ['',' ',None,'INDEF'] and hdrext is not None:
        selections = {'hdrname':fobj[hdrext].header['hdrname']}
    else:
        selections = {'distname':distname}
    wcscorr.delete_wcscorr_row(fobj['WCSCORR'].data,selections)

    # delete the headerlet extension now
    for hdrind in hdrlet_ind:
        del fobj[hdrind]
    
    # Update file object with changes
    fobj.flush()
    # close file, if was opened by this function
    if close_fobj:
        fobj.close()
    print 'Deleted headerlet from extension(s): ',hdrlet_ind

def headerlet_summary(filename,columns=None,pad=2,maxwidth=None, 
                        output=None,clobber=True,quiet=False):
    """
    Print a summary of all HeaderletHDUs in a science file to STDOUT, and
    optionally to a text file   
    The summary includes:   
        HDRLET_ext_number  HDRNAME  WCSNAME DISTNAME SIPNAME NPOLFILE D2IMFILE

    Parameters
    ----------
    filename: string or HDUList
           Either a filename or PyFITS HDUList object for the input science file
            An input filename (str) will be expanded as necessary to interpret 
            any environmental variables included in the filename.
    columns: list
        List of headerlet PRIMARY header keywords to report in summary
        By default (set to None), it will use the default set of keywords
        defined as the global list DEFAULT_SUMMARY_COLS
    pad: int 
        Number of padding spaces to put between printed columns 
        [Default: 2]
    maxwidth: int
        Maximum column width(not counting padding) for any column in summary
        By default (set to None), each column's full width will be used
    output: string (optional)
        Name of optional output file to record summary. This filename
        can contain environment variables. 
        [Default: None]
    clobber: bool
        If True, will overwrite any previous output file of same name
    quiet: bool
        If True, will NOT report info to STDOUT
    """    
    if columns is None:
        summary_cols = DEFAULT_SUMMARY_COLS
    else:
        summary_cols = columns

    summary_dict = {}
    for kw in summary_cols:
        summary_dict[kw] = copy.deepcopy(COLUMN_DICT)

    # Define Extension number column 
    extnums_col = copy.deepcopy(COLUMN_DICT)
    extnums_col['name'] = 'EXTN'
    extnums_col['width'] = 6
        
    fobj,fname,close_fobj = parseFilename(filename)
    # find all HDRLET extensions and combine info into a single summary
    for extn in fobj:
        if 'extname' in extn.header and extn.header['extname'] == 'HDRLET':
            hdrlet_indx = fobj.index_of(('hdrlet',extn.header['extver']))
            try:
                ext_cols, ext_summary = extn.headerlet.summary(columns=summary_cols)
                extnums_col['vals'].append(hdrlet_indx)
                for kw in summary_cols:
                    for key in COLUMN_DICT:
                        summary_dict[kw][key].extend(ext_summary[kw][key])
            except:
                print 'Skipping headerlet...Could not read Headerlet from extension ',hdrlet_indx

    if close_fobj:
        fobj.close()

    # Print out the summary dictionary
    print_summary(summary_cols, summary_dict, pad=pad, maxwidth=maxwidth, 
                    idcol=extnums_col, output=output, 
                    clobber=clobber, quiet=quiet)
        
def restore_from_headerlet(filename, hdrname=None, hdrext=None, archive=True, force=False): 
    """
    Restores a headerlet as a primary WCS
    
    Parameters
    ----------
    filename: string or HDUList
           Either a filename or PyFITS HDUList object for the input science file
            An input filename (str) will be expanded as necessary to interpret 
            any environmental variables included in the filename.
    hdrname: string
        HDRNAME keyword of HeaderletHDU 
    hdrext: int or tuple
        Headerlet extension number of tuple ('HDRLET',2) 
    archive: boolean (default: True) 
        When the distortion model in the headerlet is the same as the distortion model of
        the science file, this flag indicates if the primary WCS should be saved as an alternate 
        nd a headerlet extension.
        When the distortion models do not match this flag indicates if the current primary and 
        alternate WCSs should be archived as headerlet extensions and alternate WCS.
    force: boolean (default:False)
        When the distortion models of the headerlet and the primary do not match, and archive 
        is False, this flag forces an update of the primary.
    """
    initLogging('restore_from_headerlet')

    hdrlet_ind = findHeaderletHDUs(filename,hdrext=hdrext,hdrname=hdrname)

    fobj,fname,close_fobj = parseFilename(filename,mode='update')
        
    if len(hdrlet_ind) > 1:
        if hdrext:
            kwerr = 'hdrext'
            kwval = hdrext
        else:
            kwerr = 'hdrname'
            kwval = hdrname
        print '====================================================='
        print '[restore_from_headerlet]'
        print 'Multiple Headerlet extensions found with the same name.'
        print '  %d Headerlets with "%s" = %s found in %s.'%(
                    len(hdrlet_ind),kwerr,kwval,fname)
        print '====================================================='
        if close_fobj: 
            fobj.close()
        raise ValueError
        
    hdrlet_indx = hdrlet_ind[0]
    
    # read headerlet from HeaderletHDU into memory
    if hasattr(fobj[hdrlet_ind[0]], 'hdulist'):
        hdrlet = fobj[hdrlet_indx].hdulist
    else:
        hdrlet = fobj[hdrlet_indx].headerlet # older convention in PyFITS

    # read in the names of the extensions which HeaderletHDU updates
    extlist = []
    for ext in hdrlet:
        if 'extname' in ext.header and ext.header['extname'] == 'SIPWCS':
            # convert from string to tuple or int
            sciext = eval(ext.header['sciext']) 
            extlist.append(fobj[sciext])
    # determine whether distortion is the same 
    current_distname = hdrlet[0].header['distname']
    same_dist = True
    if current_distname != fobj[0].header['distname']:
        same_dist = False 
        if not archive and not force:
            if close_fobj:
                fobj.close()
            print '====================================================='
            print '[restore_from_headerlet]'
            print 'Headerlet does not have the same distortion as image!'
            print '  Set "archive"=True to save old distortion model, or'
            print '  set "force"=True to overwrite old model with new.'
            print '====================================================='
            raise ValueError

    # check whether primary WCS has been archived already 
    # Use information from first 'SCI' extension 
    priwcs_name = None

    scihdr = extlist[0].header
    sci_wcsnames = altwcs.wcsnames(scihdr).values()      
    if 'hdrname' in scihdr:
        priwcs_hdrname = scihdr['hdrname']
    else:
        if 'wcsname' in scihdr:
            priwcs_hdrname = priwcs_name = scihdr['wcsname']
        else:
            if 'idctab' in scihdr:
                priwcs_hdrname = ''.join(['IDC_',
                        utils.extract_rootname(scihdr['idctab'])])
            else:
                priwcs_hdrname = 'UNKNOWN'
            priwcs_name = priwcs_hdrname
            scihdr.update('WCSNAME',priwcs_name)

    priwcs_unique = verifyHdrnameIsUnique(fobj,priwcs_hdrname)  
    if archive and priwcs_unique:
        if priwcs_unique:
            newhdrlet = create_headerlet(fobj,sciext=scihdr['extname'],
                        hdrname=priwcs_hdrname)
            newhdrlet.attach_to_file(fobj)
    #
    # copy hdrlet as a primary
    #
    hdrlet.apply_as_primary(fobj,attach=False, archive=archive,force=force)

    fobj.flush()    
    if close_fobj:
        fobj.close()
        
def restore_all_with_distname(filename, distname, primary, archive=True, sciext='SCI'): 
    """   
    Restores all HeaderletHDUs with a given distortion model as alternate WCSs and a primary   
    
    Parameters
    --------------
    filename: string or HDUList
           Either a filename or PyFITS HDUList object for the input science file
            An input filename (str) will be expanded as necessary to interpret 
            any environmental variables included in the filename.
    distname: string 
        distortion model as represented by a DISTNAME keyword         
    primary: int or string or None
        HeaderletHDU to be restored as primary
        if int - a fits extension
        if string - HDRNAME
        if None - use first HeaderletHDU
    archive: boolean (default True)
        flag indicating if HeaderletHDUs should be created from the
        primary and alternate WCSs in fname before restoring all matching
        headerlet extensions
    """
    initLogging('restore_all_with_distname')

    fobj,fname,close_fobj = parseFilename(filename,mode='update')
        
    hdrlet_ind = findHeaderletHDUs(fobj,distname=distname)
    if len(hdrlet_ind) == 0:
        print '====================================================='
        print '[restore_all_with_distname]'
        print 'No Headerlet extensions found with '
        print '  "DISTNAME" = %s in %s.'%(kwval,fname)
        print 'Full list of DISTNAMEs found in all headerlet extensions: '
        print getHeaderletKwNames(fobj,kw='DISTNAME')
        print '====================================================='
        if close_fobj:
            fobj.close()
        raise ValueError
    
    # Interpret 'primary' parameter input into extension number
    if primary is None:
        primary_ind = hdrlet_ind[0]
    elif isinstance(primary, int):
        primary_ind = primary
    else:
        primary_ind = None
        for ind in hdrlet_ind:
            if fobj[ind].header['hdrname'] == primary:
                primary_ind = ind
                break
        if primary_ind is None:
            if close_fobj:
                fobj.close()
            print '====================================================='
            print '[restore_all_from_distname]'
            print 'No Headerlet extensions found with '
            print '  "DISTNAME" = %s in %s.'%(primary,fname)
            print '====================================================='
            raise ValueError
    # Check to see whether 'primary' HeaderletHDU has same distname as user
    # specified on input
                
    # read headerlet from HeaderletHDU into memory
    if hasattr(fobj[primary_ind], 'hdulist'):
        primary_hdrlet = fobj[primary_ind].hdulist
    else:
        primary_hdrlet = fobj[primary_ind].headerlet # older convention in PyFITS
    pri_distname = primary_hdrlet[0].header['distname']
    if pri_distname != distname:
            if close_fobj:
                fobj.close()
            print '====================================================='
            print '[restore_all_from_distname]'
            print 'Headerlet extension to be used as PRIMARY WCS '
            print '  has "DISTNAME" = %s while  %s.'%(pri_distname)
            print '  "DISTNAME" = %s was specified on input.'%(distname)
            print '  All updated WCSs must have same DISTNAME. Quitting...'
            print '====================================================='
            raise ValueError
        
    # read in the names of the WCSs which the HeaderletHDUs will update
    wnames = altwcs.wcsnames(fobj[sciext,1].header)
    
    # work out how many HeaderletHDUs will be used to update the WCSs
    numhlt = len(hdrlet_ind)
    hdrnames = getHeaderletKwNames(fobj,kw='wcsname')
    
    # read in headerletHDUs and update WCS keywords
    for hlet in hdrlet_ind:
        if fobj[hlet].header['distname'] == distname:
            if hasattr(fobj[hlet], 'hdulist'):
                hdrlet = fobj[hlet].hdulist
            else:
                hdrlet = fobj[hlet].headerlet # older convention in PyFITS
            if hlet == primary_ind:
                hdrlet.apply_as_primary(fobj,attach=False,archive=archive,force=True)
            else:
                hdrlet.apply_as_alternate(fobj,attach=False,wcsname=hdrlet[0].header['wcsname'])

    fobj.flush()
    if close_fobj:
        fobj.close()

def archive_as_headerlet(filename, hdrname, sciext='SCI', 
                        wcsname=None, wcskey=None, destim=None,
                        sipname=None, npolfile=None, d2imfile=None, 
                        author=None, descrip=None, history=None,
                        rms_ra=None, rms_dec=None, nmatch=None, catalog=None):
    """
    Save a WCS as a headerlet extension and write it out to a file.

    This function will create a headerlet, attach it as an extension to the
    science image (if it has not already been archived) then, optionally, 
    write out the headerlet to a separate headerlet file. 
    
    Either wcsname or wcskey must be provided, if both are given, they must match a valid WCS
    Updates wcscorr if necessary.

    Parameters
    ----------
    filename: string or HDUList
           Either a filename or PyFITS HDUList object for the input science file
            An input filename (str) will be expanded as necessary to interpret 
            any environmental variables included in the filename.
    hdrname: string
        Unique name for this headerlet, stored as HDRNAME keyword 
    sciext: string
        name (EXTNAME) of extension that contains WCS to be saved
    wcsname: string
        name of WCS to be archived, if " ": stop
    wcskey: one of A...Z or " " or "PRIMARY" 
        if " " or "PRIMARY" - archive the primary WCS
    destim: string
        DESTIM keyword
        if  NOne, use ROOTNAME or science file name
    sipname: string or None (default)
             Name of unique file where the polynomial distortion coefficients were
             read from. If None, the behavior is:
             The code looks for a keyword 'SIPNAME' in the science header
             If not found, for HST it defaults to 'IDCTAB'
             If there is no SIP model the value is 'NOMODEL'
             If there is a SIP model but no SIPNAME, it is set to 'UNKNOWN'
    npolfile: string or None (default)
             Name of a unique file where the non-polynomial distortion was stored.
             If None:
             The code looks for 'NPOLFILE' in science header.
             If 'NPOLFILE' was not found and there is no npol model, it is set to 'NOMODEL'
             If npol model exists, it is set to 'UNKNOWN'
    d2imfile: string
             Name of a unique file where the detector to image correction was
             stored. If None:
             The code looks for 'D2IMFILE' in the science header.
             If 'D2IMFILE' is not found and there is no d2im correction,
             it is set to 'NOMODEL'
             If d2im correction exists, but 'D2IMFILE' is missing from science
             header, it is set to 'UNKNOWN'
    author: string
            Name of user who created the headerlet, added as 'AUTHOR' keyword
            to headerlet PRIMARY header
    descrip: string
            Short description of the solution provided by the headerlet
            This description will be added as the single 'DESCRIP' keyword
            to the headerlet PRIMARY header 
    history: filename, string or list of strings
            Long (possibly multi-line) description of the solution provided
            by the headerlet. These comments will be added as 'HISTORY' cards
            to the headerlet PRIMARY header
            If filename is specified, it will format and attach all text from
            that file as the history.
    """
    initLogging('archive_as_headerlet')

    if wcsname in [None,' ','','INDEF'] and wcskey is None:
        print '='*60
        print '[archive_as_headerlet]'
        print 'No valid WCS found found in %s.'%fname
        print '   A valid value for either "wcsname" or "wcskey" '
        print '   needs to be specified. '
        print '='*60
        raise ValueError
    
    if hdrname in [None, ' ','']:
        print '='*60
        print '[archive_as_headerlet]'
        print 'No valid name for this headerlet was provided for %s.'%fname
        print '   A valid value for "hdrname" needs to be specified. '
        print '='*60
        raise ValueError
        
    # Translate 'wcskey' value for PRIMARY WCS to valid altwcs value of ' '
    if wcskey == 'PRIMARY': wcskey = ' '
    
    fobj,fname,close_fobj = parseFilename(filename,mode='update')

    numhlt = countExtn(fobj, 'HDRLET')

    if wcsname is None:
        scihdr = fobj[sciext,1].header
        wcsname = scihdr['wcsname'+wcskey]

            
    # Check to see whether or not a HeaderletHDU with this hdrname already 
    # exists
    hdrnames = getHeaderletKwNames(fobj)
    if hdrname not in hdrnames:
        hdrletobj = getHeaderletObj(fobj,sciext=sciext,
                                    wcsname=wcsname, wcskey=wcskey,
                                    hdrname=hdrname,
                                    sipname=sipname, npolfile=npolfile, 
                                    d2imfile=d2imfile, author=author, 
                                    descrip=descrip, history=history,
                                    rms_ra=rms_ra, rms_dec=rms_dec, nmatch=nmatch, 
                                    catalog=catalog,                                    
                                    hdrletnum=numhlt + 1, verbose=False) 
        hlt_hdu = HeaderletHDU.fromheaderlet(hdrletobj)

        if destim is not None:
            hdrlet_hdu[0].header['destim'] = destim

        fobj.append(hdrlet_hdu)
        
        fobj.flush()
    else:
        print 'WARNING:'
        print '    Headerlet with hdrname ',hdrname,' already archived for WCS ',wcsname
        print '    No new headerlet appended to ',fname,'.'
        
    if close_fobj:
        fobj.close()
    
#### Headerlet Class definitions
class Headerlet(pyfits.HDUList):
    """
    A Headerlet class
    Ref: http://mediawiki.stsci.edu/mediawiki/index.php/Telescopedia:Headerlets
    """

    def __init__(self, fobj, mode='copyonwrite', verbose=False, logmode='w'):
        """
        Parameters
        ----------
        fobj:  string
                Name of headerlet file, file-like object, a list of HDU
                instances, or an HDUList instance
        mode: string, optional
                Mode with which to open the given file object
        verbose: int 
                python logging level, higher numbers trigger more output
        logmode: 'w' or 'a'
                for internal use only, indicates whether the log file 
                should be open in attach or write mode
        """
        self.verbose = verbose
        self.hdr_logger = logging.getLogger('headerlet.Headerlet')
        initLogging('class Headerlet', logger=self.hdr_logger, level=100, 
                    verbose=self.verbose, logmode=logmode)

        fobj,fname,close_file = parseFilename(fobj)
        
        super(Headerlet, self).__init__(fobj)
        self.fname = self.filename()
        self.hdrname = self[0].header["HDRNAME"]
        self.wcsname = self[0].header["WCSNAME"]
        self.upwcsver = self[0].header.get("UPWCSVER", "")
        self.pywcsver = self[0].header.get("PYWCSVER", "")
        self.destim = self[0].header["DESTIM"]
        self.sipname = self[0].header["SIPNAME"]
        self.npolfile = self[0].header["NPOLFILE"]
        self.d2imfile = self[0].header["D2IMFILE"]
        self.distname = self[0].header["DISTNAME"]
        self.vafactor = self[1].header.get("VAFACTOR", 1) #None instead of 1?
        self.author = self[0].header["AUTHOR"]
        self.descrip = self[0].header["DESCRIP"]

        self.history = ''
        for card in self[0].header['HISTORY*']:
            self.history += card.value+'\n'
        
        self.d2imerr = 0
        self.axiscorr = 1

    def apply_as_primary(self, fobj, attach=True, archive=True, force=False):
        """
        Copy this headerlet as a primary WCS to fobj
        
        Parameters
        ----------
        fobj: string, HDUList
              science file to which the headerlet should be applied
        attach: boolean
              flag indicating if the headerlet should be attached as a 
              HeaderletHDU to fobj. If True checks that HDRNAME is unique 
              in the fobj and stops if not.
        archive: boolean (default is True)
              When the distortion model in the headerlet is the same as the 
              distortion model of the science file, this flag indicates if 
              the primary WCS should be saved as an alternate and a headerlet 
              extension.
              When the distortion models do not match this flag indicates if 
              the current primary and alternate WCSs should be archived as 
              headerlet extensions and alternate WCS.
        force: boolean (default is False)
              When the distortion models of the headerlet and the primary do 
              not match, and archive is False this flag forces an update 
              of the primary
        """
        self.hverify()
        if self.verify_dest(fobj):
            if not isinstance(fobj, pyfits.HDUList):
                fobj = pyfits.open(fobj, mode='update')
                close_dest = True
            else:
                close_dest = False

            # Check to see whether the distortion model in the destination
            # matches the distortion model in the headerlet being applied
            dist_models_equal=True
            if  self[0].header['DISTNAME'] != fobj[0].header['DISTNAME']:
                if self.verbose:
                    print 'Distortion model in headerlet not the same as destination model'
                    print '    Headerlet model  : ',self[0].header['DISTNAME']
                    print '    Destination model: ',fobj[0].header['DISTNAME']
                dist_models_equal = False

            if not dist_models_equal and not force:
                raise ValueError
            
            orig_hlt_hdu = None
            numhlt = countExtn(fobj, 'HDRLET')
            hdrlet_extnames = getHeaderletKwNames(fobj)
            
            # Insure that WCSCORR table has been created with all original
            # WCS's recorded prior to adding the headerlet WCS
            wcscorr.init_wcscorr(fobj)

            alt_hlethdu = []
            # If archive has been specified
            #   regardless of whether or not the distortion models are equal...
            if archive:
                if 'wcsname' in fobj[('SCI',1)].header: 
                    hdrname = fobj[('SCI',1)].header['WCSNAME']
                    wcsname = hdrname
                else:
                    hdrname = fobj[0].header['ROOTNAME'] + '_orig'
                    wcsname = None
                wcskey = ' '
                # Check the HDRNAME for all current headerlet extensions
                # to see whether this PRIMARY WCS has already been appended
                if hdrname not in hdrlet_extnames:
                    # -  if WCS has not been saved, write out WCS as headerlet extension
                    # Create a headerlet for the original Primary WCS data in the file,
                    # create an HDU from the original headerlet, and append it to
                    # the file
                    orig_hlt = getHeaderletObj(fobj,sciext='SCI',
                                    wcsname=wcsname, wcskey=wcskey,
                                    hdrname=hdrname, sipname=None, 
                                    npolfile=None, d2imfile=None, 
                                    author=None, descrip=None, history=None,
                                    hdrletnum=numhlt + 1, verbose=self.verbose) 
                    orig_hlt_hdu = HeaderletHDU.fromheaderlet(orig_hlt)
                    numhlt += 1
                    orig_hlt_hdu.header.update('EXTVER',numhlt)

                wcsextn = mapFitsExt2HDUListInd(fobj.filename(),"SCI")[('SCI',1)]
                if dist_models_equal:
                    # Use the WCSNAME to determine whether or not to archive 
                    # Primary WCS as altwcs
                    # wcsname = hwcs.wcs.name 
                    scihdr = fobj[wcsextn].header
                    if 'hdrname' in scihdr:
                        priwcs_name = scihdr['hdrname']
                    else:
                        if 'wcsname' in scihdr:
                            priwcs_name = scihdr['wcsname']
                        else:
                            if 'idctab' in scihdr:
                                priwcs_name = ''.join(['IDC_',
                                        utils.extract_rootname(scihdr['idctab'])])
                            else:
                                priwcs_name = 'UNKNOWN'
                    nextkey = altwcs.next_wcskey(fobj,ext=wcsextn)
                    numsci = countExtn(fobj,'SCI')
                    sciext_list = []
                    for i in range(1,numsci+1):
                        sciext_list.append(('SCI',i))
                    altwcs.archiveWCS(fobj,ext=sciext_list,wcskey=nextkey,wcsname=priwcs_name)
                else:
                    for hname in altwcs.wcsnames(fobj,ext=wcsextn).values():
                        if hname != 'OPUS' and hname not in hdrlet_extnames:
                            # get HeaderletHDU for alternate WCS as well
                            alt_hlet = getHeaderletObj(fobj, sciext='SCI',
                                    wcsname=hname, wcskey=wcskey,
                                    hdrname=hname, sipname=None, 
                                    npolfile=None, d2imfile=None, 
                                    author=None, descrip=None, history=None,
                                    hdrletnum=numhlt + 1, verbose=self.verbose) 
                            numhlt += 1
                            alt_hlet_hdu = HeaderletHDU.fromheaderlet(alt_hlet)
                            alt_hlet_hdu.header.update('EXTVER',numhlt)
                            alt_hlethdu.append(alt_hlet_hdu)
                            hdrlet_extnames.append(hname)

            if not dist_models_equal:
                self._delDestWCS(fobj)
                #! Always attach these extensions last. 
                # Otherwise their headers may get updated with the other WCS kw.
                numwdvar = countExtn(self, 'WCSDVARR')
                numd2im = countExtn(self, 'D2IMARR')
                for idx in range(1, numwdvar + 1):
                    fobj.append(self[('WCSDVARR', idx)].copy())
                for idx in range(1, numd2im + 1):
                    fobj.append(self[('D2IMARR', idx)].copy())

            refs = updateRefFiles(self[0].header.ascard, fobj[0].header.ascard, verbose=self.verbose)
            numsip = countExtn(self, 'SIPWCS')
            for idx in range(1, numsip + 1):
                fhdr = fobj[('SCI', idx)].header
                siphdr = self[('SIPWCS', idx)].header.ascard                

                if dist_models_equal:
                    hwcs = HSTWCS(fobj,ext=('SCI',idx))
                    hwcshdr = hwcs.wcs2header(sip2hdr=not(dist_models_equal))

                # a minimal attempt to get the position of the WCS keywords group
                # in the header by looking for the PA_APER kw.
                # at least make sure the WCS kw are written before the HISTORY kw
                # if everything fails, append the kw to the header
                akeywd = None
                bkeywd = None
                if 'PA_APER' in fhdr:
                    akeywd = 'PA_APER'
                else:
                    if 'HISTORY' in fhdr:
                        bkeywd = 'HISTORY'
                self.hdr_logger.debug(
                    "Updating WCS keywords after %s and/or before %s " %
                    (akeywd,bkeywd))
                update_cpdis = False
                for k in siphdr[-1::-1]:
                    # Replace or add WCS keyword from headerlet as PRIMARY WCS
                    # In the case that the distortion models are not equal,
                    # this will copy all keywords from headerlet into fobj
                    # When the distortion models are equal, though, it will
                    # only copy the primary WCS keywords (CRVAL,CRPIX,...) 
                    if (dist_models_equal and (k.key in hwcshdr)) or \
                     (not dist_models_equal and k.key not in FITS_STD_KW):
                        if 'DP' not in k.key:
                            fhdr.update(k.key,k.value,comment=k.comment,
                                        after=akeywd,before=bkeywd)
                        else:
                            update_cpdis = True
                    else:
                        pass
                # Update WCS with HDRNAME as well
                fhdr.update('HDRNAME',self[0].header['hdrname'],after='WCSNAME')
                
                # Update header with record-valued keywords here
                if update_cpdis:
                    numdp = len(siphdr['CPDIS*'])
                    for dpaxis in range(1,numdp+1):
                        cpdis_indx = fhdr.ascard.index_of('CPDIS%d'%(dpaxis))
                        for dpcard in siphdr['DP%d*'%(dpaxis)][-1::-1]:
                            fhdr.ascard.insert(cpdis_indx,dpcard)

            # Update the WCSCORR table with new rows from the headerlet's WCSs
            wcscorr.update_wcscorr(fobj, self, 'SIPWCS')

            # Append the original headerlet
            if archive and orig_hlt_hdu:
                fobj.append(orig_hlt_hdu)
            # Append any alternate WCS Headerlets
            if len(alt_hlethdu) > 0:
                for ahdu in alt_hlethdu:
                    fobj.append(ahdu)
            if attach:
                # Finally, append an HDU for this headerlet
                new_hlt = HeaderletHDU.fromheaderlet(self)
                new_hlt.update_ext_version(numhlt + 1)
                fobj.append(new_hlt)

            if close_dest:
                fobj.close()
        else:
            self.hdr_logger.critical("Observation %s cannot be updated with headerlet "
                            "%s" % (fobj.filename(), self.hdrname))
            print "Observation %s cannot be updated with headerlet %s" \
                  % (fobj.filename(), self.hdrname)

    def apply_as_alternate(self, fobj, attach=True, wcskey=None, wcsname=None):
        """
        Copy this headerlet as an alternate WCS to fobj
        
        Parameters
        ----------
        fobj: string, HDUList
              science file/HDUList to which the headerlet should be applied
        attach: boolean
              flag indicating if the headerlet should be attached as a 
              HeaderletHDU to fobj. If True checks that HDRNAME is unique 
              in the fobj and stops if not.
        wcskey: string
              Key value (A-Z, except O) for this alternate WCS
              If None, the next available key will be used 
        wcsname: string
              Name to be assigned to this alternate WCS
              WCSNAME is a required keyword in a Headerlet but this allows the
              user to change it as desired.
        
                    """
        self.hverify()
        if self.verify_dest(fobj):
            if not isinstance(fobj, pyfits.HDUList):
                fobj = pyfits.open(fobj, mode='update')
                close_dest = True
            else:
                close_dest = False
            fname = fobj.filename()
            
            # Verify whether this headerlet has the same distortion found in
            # the image being updated
            if 'DISTNAME' in fobj[0].header:
                distname = fobj[0].header['DISTNAME']
            else:
                # perhaps call 'updatewcs.utils.construct_distname()' instead
                distname = 'UNKNOWN' 
                
            if distname == 'UNKNOWN' or self.distname != distname:
                self.hdr_logger.critical("Observation %s cannot be updated with headerlet %s"\
                                        % (fname, self.hdrname))
                print ("Observation %s cannot be updated with headerlet %s"\
                                % (fname, self.hdrname))
                self.hdr_logger.critical("  Distortion in image:  %s \n    did not match \n headerlet distortion: %s"
                                % (distname, self.distname))
                print ( "  Distortion in image:  %s \n  did not match \n headerlet distortion: %s"
                                % (distname, self.distname))
                print ("The method .attach_to_file() can be used to append this headerlet to %s"\
                                % (fname))
                if close_dest: fobj.close()
                raise ValueError

            # Insure that WCSCORR table has been created with all original
            # WCS's recorded prior to adding the headerlet WCS
            wcscorr.init_wcscorr(fobj)
 
            # determine value of WCSNAME to be used
            if wcsname is not None:
                wname = wcsname
            else:
                wname = self[0].header['WCSNAME']

            numhlt = countExtn(fobj, 'HDRLET')
            numsip = countExtn(self, 'SIPWCS')            
            for idx in range(1, numsip + 1):
                fhdr = fobj[('SCI', idx)].header
                siphdr = self[('SIPWCS', idx)].header.ascard
                
                # determine what alternate WCS this headerlet will be assigned to
                if wcskey is None:
                    wkey = altwcs.next_wcskey(fobj[('SCI',idx)].header)
                else:
                    available_keys = altwcs.available_wcskeys(fobj[('SCI',idx)].header)
                    if wcskey in available_keys:
                        wkey = wcskey
                    else:
                        self.hdr_logger.critical("Observation %s already contains alternate WCS with key \
                                        %s" % (fname, wcskey))
                        print ("Observation %s already contains alternate WCS with key \
                                        %s" % (fname, wcskey))
                        if close_dest: fobj.close()
                        raise ValueError

                # a minimal attempt to get the position of the WCS keywords group
                # in the header by looking for the PA_APER kw.
                # at least make sure the WCS kw are written before the HISTORY kw
                # if everything fails, append the kw to the header
                try:
                    wind = fhdr.ascard.index_of('HISTORY')
                except KeyError:
                    wind = len(fhdr)
                self.hdr_logger.debug("Inserting WCS keywords at index %s" % wind)

                for k in siphdr:
                    """
                    if k.key not in ['XTENSION', 'BITPIX', 'NAXIS', 'PCOUNT',
                                     'GCOUNT','EXTNAME', 'EXTVER', 'ORIGIN',
                                     'INHERIT', 'DATE', 'IRAF-TLM']:
                    """
                    for akw in altwcs.altwcskw:
                        if akw in k.key:
                            fhdr.ascard.insert(wind,pyfits.Card(
                                            key=k.key[:7]+wkey,value=k.value,
                                            comment=k.comment))
                    else:
                        pass

                fhdr.ascard.insert(wind,pyfits.Card('WCSNAME'+wkey,wname))
                # also update with HDRNAME (a non-WCS-standard kw)
                fhdr.ascard.insert(wind,pyfits.Card('HDRNAME'+wkey,
                                        self[0].header['hdrname']))

            # Update the WCSCORR table with new rows from the headerlet's WCSs
            wcscorr.update_wcscorr(fobj, self, 'SIPWCS')
            
            if attach:
                # Finally, append an HDU for this headerlet
                new_hlt = HeaderletHDU.fromheaderlet(self)
                new_hlt.update_ext_version(numhlt + 1)
                fobj.append(new_hlt)

            if close_dest:
                fobj.close()
        else:
            self.hdr_logger.critical("Observation %s cannot be updated with headerlet "
                            "%s" % (fname, self.hdrname))
            print "Observation %s cannot be updated with headerlet %s" \
                  % (fname, self.hdrname)

    def attach_to_file(self,fobj):
        """
        Attach Headerlet as an HeaderletHDU to a science file
        
        Parameters
        ----------
        fobj: string, HDUList
              science file/HDUList to which the headerlet should be applied

        Notes
        -----
        The algorithm used by this method:
        - verify headerlet can be applied to this file (based on DESTIM)
        - verify that HDRNAME is unique for this file
        - attach as HeaderletHDU to fobj
        - update wcscorr 
        """
        self.hverify()
        if not isinstance(fobj, pyfits.HDUList):
            fobj = pyfits.open(fobj, mode='update')
            close_dest = True
        else:
            close_dest = False
        if self.verify_dest(fobj) and self.verify_hdrname(fobj):

            numhlt = countExtn(fobj, 'HDRLET')
            new_hlt = HeaderletHDU.fromheaderlet(self)
            new_hlt.header.update('extver',numhlt + 1)
            fobj.append(new_hlt)

            wcscorr.update_wcscorr(fobj, self, 'SIPWCS',active=False)

        else:
            self.hdr_logger.critical("Observation %s cannot be updated with headerlet "
                            "%s" % (fobj.filename(), self.hdrname))
            print "Observation %s cannot be updated with headerlet %s" \
                  % (fobj.filename(), self.hdrname)

        if close_dest:
            fobj.close()

    def info(self, columns=None, pad=2, maxwidth=None, 
                output=None, clobber=True, quiet=False):
        """
        Prints a summary of this headerlet
        The summary includes:   
            HDRNAME  WCSNAME DISTNAME SIPNAME NPOLFILE D2IMFILE
        
        Parameters
        ----------
        columns: list
            List of headerlet PRIMARY header keywords to report in summary
            By default (set to None), it will use the default set of keywords
            defined as the global list DEFAULT_SUMMARY_COLS
        pad: int 
            Number of padding spaces to put between printed columns 
            [Default: 2]
        maxwidth: int
            Maximum column width(not counting padding) for any column in summary
            By default (set to None), each column's full width will be used
        output: string (optional)
            Name of optional output file to record summary. This filename
            can contain environment variables. 
            [Default: None]
        clobber: bool
            If True, will overwrite any previous output file of same name
        quiet: bool
            If True, will NOT report info to STDOUT

        """
        summary_cols, summary_dict = self.summary(columns=columns)        
        print_summary(summary_cols, summary_dict, pad=pad, maxwidth=maxwidth, 
                        idcol=None, output=output, clobber=clobber, quiet=quiet)

    def summary(self, columns=None):
        """
        Returns a summary of this headerlet as a dictionary
        
        The summary includes a summary of the distortion model as :   
            HDRNAME  WCSNAME DISTNAME SIPNAME NPOLFILE D2IMFILE
        
        Parameters
        ----------
        columns: list
            List of headerlet PRIMARY header keywords to report in summary
            By default(set to None), it will use the default set of keywords
            defined as the global list DEFAULT_SUMMARY_COLS
            
        Returns
        -------
        summary: dict
            Dictionary of values for summary
        """
        if columns is None:
            summary_cols = DEFAULT_SUMMARY_COLS
        else:
            summary_cols = columns

        # Initialize summary dict based on requested columns
        summary = {}
        for kw in summary_cols:
            summary[kw] = copy.deepcopy(COLUMN_DICT)

        # Populate the summary with headerlet values
        for kw in summary_cols:
            if kw in self[0].header:
                val = self[0].header[kw]
            else:
                val = 'INDEF'
            summary[kw]['vals'].append(val)
            summary[kw]['width'].append(max(len(val),len(kw)))
        
        return summary_cols,summary

    def hverify(self):
        self.verify()
        header = self[0].header
        assert('DESTIM' in header and header['DESTIM'].strip())
        assert('HDRNAME' in header and header['HDRNAME'].strip())
        assert('STWCSVER' in header)
       
    def verify_hdrname(self,dest):
        """
        Verifies that the headerlet can be applied to the observation
        
        Reports whether or not this file already has a headerlet with this 
        HDRNAME.
        """
        unique = verifyHdrnameIsUnique(dest,self.hdrname)
        self.hdr_logger.debug("verify_hdrname() returned %s"%unique)
        return unique
    
    def verify_dest(self, dest):
        """
        verifies that the headerlet can be applied to the observation

        DESTIM in the primary header of the headerlet must match ROOTNAME
        of the science file (or the name of the destination file)
        """

        try:
            if not isinstance(dest, pyfits.HDUList):
                droot = pyfits.getval(dest, 'ROOTNAME')
            else:
                droot = dest[0].header['ROOTNAME']
        except KeyError:
            self.hdr_logger.debug("Keyword 'ROOTNAME' not found in destination file")
            droot = dest.split('.fits')[0]
        if droot == self.destim:
            self.hdr_logger.debug("verify_destim() returned True")
            return True
        else:
            self.hdr_logger.debug("verify_destim() returned False")
            return False

    def tofile(self, fname, destim=None, hdrname=None, clobber=False):
        if not destim or not hdrname:
            self.hverify()
        self.writeto(fname, clobber=clobber)

    def _delDestWCS(self, dest):
        """
        Delete the WCS of a science file
        """

        self.hdr_logger.info("Deleting all WCSs of file %s" % dest.filename())
        numext = len(dest)

        for idx in range(numext):
            # Only delete WCS from extensions which may have WCS keywords
            if ('XTENSION' in dest[idx].header and
                dest[idx].header['XTENSION'] == 'IMAGE'):
                self._removeD2IM(dest[idx])
                self._removeSIP(dest[idx])
                self._removeLUT(dest[idx])
                self._removePrimaryWCS(dest[idx])
                self._removeIDCCoeffs(dest[idx])
                try:
                    del dest[idx].header.ascard['VAFACTOR']
                except KeyError:
                    pass

        self._removeRefFiles(dest[0])
        self._removeAltWCS(dest, ext=range(numext))
        numwdvarr = countExtn(dest, 'WCSDVARR')
        numd2im = countExtn(dest, 'D2IMARR')
        for idx in range(1, numwdvarr + 1):
            del dest[('WCSDVARR', idx)]
        for idx in range(1, numd2im + 1):
            del dest[('D2IMARR', idx)]

    def _removeRefFiles(self, phdu):
        """
        phdu: Primary HDU
        """
        refkw = ['IDCTAB', 'NPOLFILE', 'D2IMFILE']
        for kw in refkw:
            try:
                del phdu.header.ascard[kw]
            except KeyError:
                pass

    def _removeSIP(self, ext):
        """
        Remove the SIP distortion of a FITS extension
        """

        self.hdr_logger.debug("Removing SIP distortion from (%s, %s)"
                     % (ext.name, ext._extver))
        for prefix in ['A', 'B', 'AP', 'BP']:
            try:
                order = ext.header[prefix + '_ORDER']
                del ext.header[prefix + '_ORDER']
            except KeyError:
                continue
            for i in range(order + 1):
                for j in range(order + 1):
                    key = prefix + '_%d_%d' % (i, j)
                    try:
                        del ext.header[key]
                    except KeyError:
                        pass
        try:
            del ext.header['IDCTAB']
        except KeyError:
            pass

    def _removeLUT(self, ext):
        """
        Remove the Lookup Table distortion of a FITS extension
        """

        self.hdr_logger.debug("Removing LUT distortion from (%s, %s)"
                     % (ext.name, ext._extver))
        try:
            cpdis = ext.header['CPDIS*']
        except KeyError:
            return
        try:
            for c in range(1, len(cpdis) + 1):
                del ext.header['DP%s*...' % c]
                del ext.header[cpdis[c - 1].key]
            del ext.header['CPERR*']
            del ext.header['NPOLFILE']
            del ext.header['NPOLEXT']
        except KeyError:
            pass

    def _removeD2IM(self, ext):
        """
        Remove the Detector to Image correction of a FITS extension
        """

        self.hdr_logger.debug("Removing D2IM correction from (%s, %s)"
                     % (ext.name, ext._extver))
        d2imkeys = ['D2IMFILE', 'AXISCORR', 'D2IMEXT', 'D2IMERR']
        for k in d2imkeys:
            try:
                del ext.header[k]
            except KeyError:
                pass

    def _removeAltWCS(self, dest, ext):
        """
        Remove Alternate WCSs of a FITS extension.
        A WCS with wcskey 'O' is never deleted.
        """
        dkeys = altwcs.wcskeys(dest[('SCI', 1)].header)
        self.hdr_logger.debug("Removing alternate WCSs with keys %s from %s"
                     % (dkeys, dest.filename()))
        for k in dkeys:
            if k not in ['O',' ','']: # Never delete WCS with wcskey='O'
                altwcs.deleteWCS(dest, ext=ext, wcskey=k)

    def _removePrimaryWCS(self, ext):
        """
        Remove the primary WCS of a FITS extension
        """

        self.hdr_logger.debug("Removing Primary WCS from (%s, %s)"
                     % (ext.name, ext._extver))
        naxis = ext.header.ascard['NAXIS'].value
        for key in basic_wcs:
            for i in range(1, naxis + 1):
                try:
                    del ext.header.ascard[key + str(i)]
                except KeyError:
                    pass
        try:
            del ext.header.ascard['WCSAXES']
        except KeyError:
            pass

    def _removeIDCCoeffs(self, ext):
        """
        Remove IDC coefficients of a FITS extension
        """

        self.hdr_logger.debug("Removing IDC coefficient from (%s, %s)"
                     % (ext.name, ext._extver))
        coeffs = ['OCX10', 'OCX11', 'OCY10', 'OCY11', 'IDCSCALE']
        for k in coeffs:
            try:
                del ext.header.ascard[k]
            except KeyError:
                pass


class HeaderletHDU(pyfits.hdu.base.NonstandardExtHDU):
    """
    A non-standard extension HDU for encapsulating Headerlets in a file.  These
    HDUs have an extension type of HDRLET and their EXTNAME is derived from the
    Headerlet's HDRNAME.

    The data itself is a tar file containing a single file, which is the
    Headerlet file itself.  The file name is derived from the HDRNAME keyword,
    and should be in the form `<HDRNAME>_hdr.fits`.  If the COMPRESS keyword
    evaluates to `True`, the tar file is compressed with gzip compression.

    The Headerlet contained in the HDU's data can be accessed by the
    `headerlet` attribute.
    """

    _extension = 'HDRLET'

    @pyfits.util.lazyproperty
    def data(self):
        self._file.seek(self._datLoc)
        return self._file.readarray(self.size)

    @pyfits.util.lazyproperty
    def headerlet(self):
        self._file.seek(self._datLoc)
        s = StringIO()
        # Read the data into a StringIO--reading directly from the file
        # won't work (at least for gzipped files) due to problems deep
        # within the gzip module that make it difficult to read gzip files
        # embedded in another file
        s.write(self._file.read(self.size))
        s.seek(0)
        if self._header['COMPRESS']:
            mode = 'r:gz'
        else:
            mode = 'r'
        t = tarfile.open(mode=mode, fileobj=s)
        members = t.getmembers()
        if not len(members):
            raise ValueError('The Headerlet contents are missing.')
        elif len(members) > 1:
            warnings.warn('More than one file is contained in this '
                          'only the headerlet file should be present.')
        hlt_name = self._header['HDRNAME'] + '_hdr.fits'
        try:
            hlt_info = t.getmember(hlt_name)
        except KeyError:
            warnings.warn('The file %s was missing from the HDU data.  '
                          'Assuming that the first file in the data is '
                          'headerlet file.' % hlt_name)
            hlt_info = members[0]

        hlt_file = t.extractfile(hlt_info)
        # hlt_file is a file-like object
        return Headerlet(hlt_file, mode='readonly')
    
    @classmethod
    def fromheaderlet(cls, headerlet, compress=False):
        """
        Creates a new HeaderletHDU from a given Headerlet object.

        Parameters
        ----------
        headerlet : Headerlet
            A valid Headerlet object.
        compress : bool, optional
            Gzip compress the headerlet data.
        """

        phdu = headerlet[0]
        phduhdr = phdu.header
        hlt_filename = phdu.header['HDRNAME'] + '_hdr.fits'

        # TODO: As it stands there's no good way to write out an HDUList in
        # memory, since it automatically closes the given file-like object when
        # it's done writing.  I'd argue that if passed an open file handler it
        # should not close it, but for now we'll have to write to a temp file.
        fd, name = tempfile.mkstemp()
        try:
            f = os.fdopen(fd, 'rb+')
            headerlet.writeto(f)
            # The tar file itself we'll write in memory, as it should be
            # relatively small
            if compress:
                mode = 'w:gz'
            else:
                mode = 'w'
            s = StringIO()
            t = tarfile.open(mode=mode, fileobj=s)
            t.add(name, arcname=hlt_filename)
            t.close()
        finally:
            os.remove(name)

        if 'sipname' in headerlet[0].header:
            sipname = headerlet[0].header['sipname']
        else:
            sipname = headerlet[0].header['wcsname']

        cards = [
            pyfits.Card('XTENSION', cls._extension, 'Headerlet extension'),
            pyfits.Card('BITPIX',    8, 'array data type'),
            pyfits.Card('NAXIS',     1, 'number of array dimensions'),
            pyfits.Card('NAXIS1',    len(s.getvalue()), 'Axis length'),
            pyfits.Card('PCOUNT',    0, 'number of parameters'),
            pyfits.Card('GCOUNT',    1, 'number of groups'),
            pyfits.Card('EXTNAME', cls._extension,
                        'name of the headerlet extension'),
            phdu.header.ascard['HDRNAME'],
            phdu.header.ascard['DATE'],
            pyfits.Card('SIPNAME', sipname,
                        'SIP distortion model name'),
            pyfits.Card('WCSNAME', headerlet[0].header['WCSNAME'],
                        'WCS name'),
            pyfits.Card('DISTNAME', headerlet[0].header['DISTNAME'],
                        'Distortion model name'),
            phdu.header.ascard['NPOLFILE'],
            phdu.header.ascard['D2IMFILE'],
            pyfits.Card('COMPRESS', compress, 'Uses gzip compression')
        ]

        header = pyfits.Header(pyfits.CardList(cards))

        hdu = cls(data=pyfits.DELAYED, header=header)
        # TODO: These hacks are a necessary evil, but should be removed once
        # the pyfits API is improved to better support specifying a backing
        # file object for an HDU
        hdu._file = pyfits.file._File(s)
        hdu._datLoc = 0
        return hdu

    @classmethod
    def match_header(cls, header):
        """
        This is a class method used in the pyfits refactoring branch to
        recognize whether or not this class should be used for instantiating
        an HDU object based on values in the header.

        It is included here for forward-compatibility.
        """

        card = header.ascard[0]
        if card.key != 'XTENSION':
            return False
        xtension = card.value.rstrip()
        return xtension == cls._extension

    # TODO: Add header verification

    def _summary(self):
        # TODO: Perhaps make this more descriptive...
        return (self.name, self.__class__.__name__, len(self._header.ascard))

pyfits.register_hdu(HeaderletHDU)
