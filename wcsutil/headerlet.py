from __future__ import division
import time
import numpy as np
import pyfits
from hstwcs import HSTWCS
import altwcs
from pyfits import HDUList
from mappings import basic_wcs

import logging, time
logger = logging.getLogger('stwcs.wcsutil.headerlet')

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
    if verbose == False:
        logger.setLevel(100)
    else:
        logger.setLevel(verbose)
    logger.info("Starting isWCSIdentical: %s" %time.asctime())
    
    result = True
    numsci1 = max(countext(scifile, 'SCI'), countext(scifile, 'SIPWCS'))
    numsci2 = max(countext(file2, 'SCI'), countext(file2, 'SIPWCS'))
    
    if numsci1 == 0 or numsci2==0 or numsci1!= numsci2:
        logger.info("Number of SCI and SIPWCS extensions do not match.")
        result = False
    
    if getRootname(scifile) != getRootname(file2):
        logger.info('Rootnames do not match.')
        result = False
    try:
        extname1 = pyfits.getval(scifile, 'EXTNAME', ext=('SCI', 1))
    except KeyError:
        extname1 = 'SIPWCS'
    try:
        extname2 = pyfits.getval(file2, 'EXTNAME', ext=('SCI', 1))
    except KeyError:
        extname2 = 'SIPWCS'    
        
    for i in range(1, numsci1+1):
        w1 = HSTWCS(scifile,ext=(extname1,i))
        w2 = HSTWCS(file2, ext=(extname2,i))
        if not (w1.wcs.crval == w2.wcs.crval).all() or \
            not (w1.wcs.crpix == w2.wcs.crpix).all()  or \
            not  (w1.wcs.cd == w2.wcs.cd).all() or \
            not (np.array(w1.wcs.ctype) == np.array(w2.wcs.ctype)).all():
            logger.info('Primary WCSs do not match')
            result = False
        if w1.sip or w2.sip: 
            if (w2.sip and not w1.sip) or (w1.sip and not w2.sip) or \
                not (w1.sip.a == w2.sip.a).all() or \
                not (w1.sip.b == w2.sip.b).all():
                logger.info('SIP coefficients do not match')
                result = False
        if w1.cpdis1 or w2.cipdis1:
            if w1.cpdis1 and not w2.cpdis1 or \
                w2.cpdis1 and not w1.cpdis1 or \
                not (w1.cpdis1.data == w2.cpdis1.data).all():
                logger.info('NPOL distortions do not match')
                result = False
        if w1.cpdis2 or w2.cpdis2:
            if w1.cpdis2 and not w2.cpdis2 or \
                w2.cpdis2 and not w1.cpdis2 or \
                not (w1.cpdis2.data == w2.cpdis2.data).all():
                logger.info('NPOL distortions do not match')
                result = False
        if w1.det2im1 or w2.det2im1:
            if w1.det2im1 and not w2.det2im1 or \
                w2.det2im1 and not w1.det2im1 or\
                not (w1.det2im1.data == w2.det2im1.data).all():
                logger.info('Det2Im corrections do not match')
                result =  False
        if w1.det2im2 or w2.det2im2:
            if w1.det2im2 and not w2.det2im2 or \
                w2.det2im2 and not w1.det2im2 or\
                not (w1.det2im2.data == w2.det2im2.data).all():
                logger.info('Det2Im corrections do not match')
                result = False
        if w1.vafactor != w2.vafactor:
            logger.info('VA factors do not match')
            result = False
    
    return result
        
def createHeaderlet(fname, hdrname, destim=None, output=None, verbose=False):
    """
    Create a headerlet from a science observation
    
    Parameters
    ----------
    fname: string
           Name of file with science observation
    hdrname: string
           Name for the headerlet, stored in the primary header of the headerlet
           If output is None, hdrname is used as output
    destim: string
           Destination image, stored in the primary header of the headerlet. 
           If None ROOTNAME is used of the science observation is used.
           ROOTNAME has precedence, destim is used for observations without 
            ROOTNAME in the primary header
    output: string
           Name for the headerlet file.
           HDRNAME is used if output is not given.
    verbose: False or a python logging level
             (one of 'INFO', 'DEBUG' logging levels) 
             (an integer representing a logging level) 
    """
    if verbose == False:
        logger.setLevel(100)
    else:
        logger.setLevel(verbose)
    logger.info("Starting createHeaderlet: %s" %time.asctime())
    fmt="%Y-%m-%dT%H:%M:%S"
    phdukw = {'IDCTAB': True, 
            'NPOLFILE': True, 
            'D2IMFILE': True} 
    fobj = pyfits.open(fname)
    if destim is None:
        try:
            destim = fobj[0].header['ROOTNAME']
        except KeyError:
            logger.exception('Required keyword "DESTIM" not found')
            print 'Please provide a value for the DESTIM keyword'
            raise
    if hdrname is None:
        logger.critical("Required keyword 'HDRNAME' not given")
        raise ValueError, "Please provide a name for the headerlet, HDRNAME is a required parameter."
    
        
        
    altkeys = altwcs.wcskeys(fobj[('SCI',1)].header)
    try:
        upwcsver = fobj[0].header.ascard['STWCSVER']
    except KeyError:
        upwcsver = pyfits.Card("STWCSVER", " ","Version of STWCS used to update the WCS")
    if 'O' in altkeys:
        altkeys.remove('O')
    numsci = countext(fname, 'SCI')
    logger.debug("Number of 'SCI' extensions in file %s is %s" % (fname, numsci))
    hdul = pyfits.HDUList()
    phdu = pyfits.PrimaryHDU()
    phdu.header.update('DESTIM', destim, comment='Destination observation root name')
    phdu.header.update('HDRNAME',hdrname, comment='Headerlet name')
    phdu.header.update('DATE',time.strftime(fmt), comment='Date FITS file was generated')
    phdu.header.ascard.append(upwcsver)
    
    updateRefFiles(fobj[0].header.ascard, phdu.header.ascard)
    phdu.header.update(key='VAFACTOR', value=fobj[('SCI',1)].header.get('VAFACTOR', 1.))
    hdul.append(phdu)
    
    for e in range(1,numsci+1):
        hwcs = HSTWCS(fname,ext=('SCI',e))
        h = hwcs.wcs2header(sip2hdr=True).ascard
        for ak in altkeys:
            awcs = HSTWCS(fname,ext=('SCI',e), wcskey=ak)
            h.extend(awcs.wcs2header(idc2hdr=False).ascard)
        h.insert(0,pyfits.Card(key='EXTNAME', value='SIPWCS', comment='Extension name'))
        h.insert(1,pyfits.Card(key='EXTVER', value=e, comment='Extension version'))
        
        fhdr = fobj[('SCI',e)].header.ascard
        if phdukw['NPOLFILE']:
            
            cpdis = fhdr['CPDIS*...']
            for c in range(1,len(cpdis)+1):
                h.append(cpdis[c-1])
                dp = fhdr['DP%s.*...'%c]
                h.extend(dp)
                
                try: h.append(fhdr['CPERROR%s'%c])
                except KeyError: pass
                
            try: h.append(fhdr['NPOLEXT'])
            except KeyError: pass
            
        if phdukw['D2IMFILE']:
            try: h.append(fhdr['D2IMEXT'])
            except KeyError: pass
            
            try: h.append(fhdr['AXISCORR'])
            except KeyError: 
                logger.exception("Required keyword AXISCORR was not found in %s['SCI',%d]" %(fname,e))
                raise
            
            try: h.append(fhdr['D2IMERR'])
            except KeyError: h.append(pyfits.Card(key='DPERROR', value=0, 
                                                comment='Maximum error of D2IMARR'))
            
        hdu = pyfits.ImageHDU(header=pyfits.Header(h))
        # temporary fix for pyfits ticket # 48
        hdu._extver = e
        hdul.append(hdu)
    numwdvarr = countext(fname, 'WCSDVARR')
    numd2im = countext(fname, 'D2IMARR')
    for w in range(1, numwdvarr+1):
        hdu = fobj[('WCSDVARR',w)].copy()
        # temporary fix for pyfits ticket # 48
        hdu._extver = w
        hdul.append(hdu)
    for d in range(1, numd2im+1):
        hdu = fobj[('D2IMARR',d)].copy()
        # temporary fix for pyfits ticket # 48
        hdu._extver = d
        hdul.append(hdu)
    if output is not None:
        # write the headerlet to a file
        if not output.endswith('_hdr.fits'):
            output = output+'_hdr.fits'
        hdul.writeto(output,clobber=True)
    fobj.close()
    return Headerlet(hdul)

def applyHeaderlet(hdrfile, destfile, destim=None, hdrname=None, createheaderlet=True, verbose=False):
    """
    Apply headerlet 'hdrfile' to a science observation 'destfile'
    
    Parameters
    ----------
    hdrfile: string
             Headerlet file
    destfile: string
             File name of science observation wyhose WCS solution willbe updated
    destim: string or None (default)
            ROOTNAME of destfile (default) 
            This string will be written as the DESTIM keyword in the headerlet, 
            created from the old WCS solution
            Normally it should be None, unless the science file is missing the 
            ROOTNAME kw, in which case the name of the file (stripped of .fits) 
            should be specified.
    hdrname: string or None (default)
            will be the value of the HDRNAME keyword in the headerlet.
            It's not  required if createheaderlet is False
    createheaderlet: boolean
            True (default): before updating, create a headerlet with the 
            WCS old solution.
    verbose: False or a python logging level
             (one of 'INFO', 'DEBUG' logging levels) 
             (an integer representing a logging level)  
    """
    logger.info("Starting applyHeaderlet: %s" %time.asctime())
    hlet = Headerlet(hdrfile)
    hlet.apply2obs(destfile, destim=destim, hdrname=hdrname, createheaderlet=createheaderlet)

def updateRefFiles(source, dest):
    """
    Update the reference files name in the primary header of 'dest'
    using values from 'source'
    
    Parameters
    ----------
    source: pyfits.Header.ascardlist
    dest:   pyfits.Header.ascardlist
    """
    logger.info("Updating reference files")
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
            phdukw[key] = False

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
    if isinstance(fname, str):
        f = pyfits.open(fname)
    else:
        f = fname
    d = {}
    for hdu in f:
        if hdu.header.has_key('EXTNAME') and hdu.header['EXTNAME'] == extname:
            extver = hdu.header['EXTVER']
            d[(extname, extver)] = f.index_of((extname, extver))
    f.close()
    return d
            
def countext(fname, extname):
    """
    Return the number of extensions with EXTNAME in the file.
    """
    if isinstance(fname, str):
        f = pyfits.open(fname)
    else:
        f = fname
    numext = 0
    for hdu in f:
        if hdu.header.has_key('EXTNAME') and hdu.header['EXTNAME'] == extname:
            numext+=1
    logger.debug("file %s has %s extensions with extname=%s" % (f.filename(), numext, extname))
    return numext

def _delDestWCS(dest):
    """
    Delete the WCS of a science file
    """
    logger.info("Deleting all WCSs of file %s" %dest) 
    fobj = pyfits.open(dest, mode='update')
    numext = len(fobj)
    
    for e in range(numext):
        _removeD2IM(fobj[e])
        _removeSIP(fobj[e])
        _removeLUT(fobj[e])
        _removePrimaryWCS(fobj[e])
        _removeIDCCoeffs(fobj[e])
        try:
            del fobj[e].header.ascard['VAFACTOR']
        except KeyError: pass
        
    _removeAltWCS(fobj, ext=range(numext)) 
    numwdvarr = countext(dest,'WCSDVARR')
    numd2im = countext(dest,'D2IMARR')
    for i in range(1, numwdvarr+1):
        del fobj[('WCSDVARR',i)]
    for i in range(1, numd2im+1):
        del fobj[('D2IMARR',i)]
    fobj.close()
    
def _removeSIP(ext):
    """
    Remove the SIP distortion of a FITS extension
    """
    logger.debug("Removing SIP distortion from (%s, %s)" % (ext.name, ext._extver))
    for prefix in ['A','B', 'AP', 'BP']:
        try:   
            order = ext.header[prefix+'_ORDER']
            del ext.header[prefix+'_ORDER']
        except KeyError:
            continue
        for i in range(order+1):
            for j in range(order+1):
                key = prefix + '_%d_%d' % (i,j)
                try:
                    del ext.header[key]
                except KeyError:
                    pass
    try: 
        del ext.header['IDCTAB']
    except KeyError:
        pass

def _removeLUT(ext):
    """
    Remove the Lookup Table distortion of a FITS extension
    """
    logger.debug("Removing LUT distortion from (%s, %s)" % (ext.name, ext._extver))
    try:
        cpdis = ext.header['CPDIS*']
    except KeyError: 
        return
    try:
        for c in range(1,len(cpdis)+1):
            del ext.header['DP%s.*...'%c]
            del ext.header[cpdis[c-1].key]   
        del ext.header['CPERR*']
        del ext.header['NPOLFILE']
        del ext.header['NPOLEXT']
    except KeyError:
        pass

def _removeD2IM(ext):
    """
    Remove the Detector to Image correction of a FITS extension
    """
    logger.debug("Removing D2IM correction from (%s, %s)" % (ext.name, ext._extver))
    d2imkeys = ['D2IMFILE', 'AXISCORR', 'D2IMEXT', 'D2IMERR']
    for k in d2imkeys:
        try:
            del ext.header[k]
        except KeyError:
            pass
        
def _removeAltWCS(dest, ext):
    """
    Remove Alternate WCSs of a FITS extension.
    A WCS with wcskey 'O' is never deleted.
    """
    dkeys = altwcs.wcskeys(dest[('SCI',1)].header)
    logger.debug("Removing alternate WCSs with keys %s from %s" % (dkeys, dest.filename()))
    for k in dkeys:
        altwcs.deleteWCS(dest, ext=ext, wcskey=k)

def _removePrimaryWCS(ext):
    """
    Remove the primary WCS of a FITS extension
    """
    logger.debug("Removing Primary WCS from (%s, %s)" % (ext.name, ext._extver))
    naxis = ext.header.ascard['NAXIS'].value
    for key in basic_wcs:
        for i in range(1,naxis+1):
            try:
                del ext.header.ascard[key+str(i)]
            except KeyError: pass
    try:
        del ext.header.ascard['WCSAXES']
    except KeyError: pass

def _removeIDCCoeffs(ext):
    """
    Remove IDC coefficients of a FITS extension
    """
    logger.debug("Removing IDC coefficient from (%s, %s)" % (ext.name, ext._extver))
    coeffs = ['OCX10', 'OCX11', 'OCY10', 'OCY11', 'IDCSCALE']
    for k in coeffs:
        try:
            del ext.header.ascard[k]
        except KeyError:
            pass
    
    
class Headerlet(HDUList):
    """
    A Headerlet class
    Ref: http://stsdas.stsci.edu/stsci_python_sphinxdocs/stwcs/headerlet_def.html
    """
    def __init__(self, fname, wcskeys=[], verbose=False):
        """
        Parameters
        ----------
        fname:  string
                Name of headerlet file
        wcskeys: python list
                a list of wcskeys to be included in the headerlet
                created from the old WCS solution before the 
                science file is updated. If empty: all alternate (if any)
                WCSs are copied to the headerlet.
        """
        logger.info("Creating a Headerlet object from wcskeys %s" % wcskeys)
        self.wcskeys = wcskeys
        if isinstance(fname,str):
            fobj = pyfits.open(fname)
        elif isinstance(fname, list):
            fobj = fname
        else:
            raise ValueError, "Input must be a file name (string) or a pyfits file object (HDUList)."
        HDUList.__init__(self, fobj) 
        self.fname = self.filename()
        self.hdrname = self[0].header["HDRNAME"]
        self.stwcsver = self[0].header.get("STWCSVER","")
        self.destim = self[0].header["DESTIM"]
        self.idctab = self[0].header.get("IDCTAB","")
        self.npolfile = self[0].header.get("NPOLFILE","")
        self.d2imfile = self[0].header.get("D2IMFILE","")
        self.vafactor = self[1].header.get("VAFACTOR", 1) #None instead of 1?
        self.d2imerr = 0
        self.axiscorr = 1
        
    def apply2obs(self,dest, destim=None, hdrname=None, createheaderlet=True):
        """
        apply this headerlet to a file
        """
        self.hverify()
        if self.verify_dest(dest):
            if createheaderlet:
                createHeaderlet(dest, hdrname, destim)
            _delDestWCS(dest)
            fobj = pyfits.open(dest, mode='update')
            
            updateRefFiles(self[0].header.ascard, fobj[0].header.ascard)
            numsip = countext(self, 'SIPWCS')
            for i in range(1,numsip+1):
                fhdr = fobj[('SCI',i)].header.ascard
                siphdr = self[('SIPWCS',i)].header.ascard
                # a minimal attempt to get the position of the WCS keywords group 
                # in the header by looking for the PA_APER kw.
                # at least make sure the WCS kw are written befir the HISTORY kw
                # if everything fails, append the kw to the header
                try:
                    wind = fhdr.index_of('PA_APER')
                except KeyError:
                    try:
                        wind = fhdr.index_of('HISTORY')
                    except KeyError:
                        wind = len(fhdr)
                logger.debug("Inserting WCS keywords at index %s" % wind)
                for k in siphdr.keys():
                    if k not in ['XTENSION', 'BITPIX', 'NAXIS', 'PCOUNT', 
                                 'GCOUNT','EXTNAME', 'EXTVER', 'ORIGIN', 
                                 'INHERIT', 'DATE', 'IRAF-TLM']:
                        fhdr.insert(wind, siphdr[k])
                    else:
                        pass
            #! Always attach these extensions last. Otherwise therr headers may
            # get updated with the other WCS kw. 
            numwdvar = countext(self,'WCSDVARR')
            numd2im = countext(self, 'D2IMARR')
            for e in range(1,numwdvar+1):
                hdu = self[('WCSDVARR',e)].copy()
                fobj.append(self[('WCSDVARR',e)].copy())
            for e in range(1, numd2im+1):
                fobj.append(self[('D2IMARR',e)].copy())
            fobj.close()
        else:
            logger.critical("Observation %s cannot be updated with headerlet %s" %(dest, self.hdrname))
            print "Observation %s cannot be updated with headerlet %s" %(dest, self.hdrname)
        
        
    def hverify(self):
        self.verify()
        assert(self[0].header.has_key('DESTIM') and self[0].header['DESTIM']!= " ")
        assert(self[0].header.has_key('HDRNAME') and self[0].header['HDRNAME']!= " ")
        assert(self[0].header.has_key('STWCSVER'))
        
    
    def verify_dest(self, dest):
        """
        verifies that the headerlet can be applied to the observation
       
        DESTIM in the primary header of the headerlet must match ROOTNAME 
        of the science file (or the name of the destination file)
        """
        try:
            droot = pyfits.getval(dest, 'ROOTNAME')
        except KeyError:
            logger.debug("Keyword 'ROOTNAME' not found in destination file")
            droot = dest.split('.fits')[0]
        if droot == self.destim:
            logger.debug("verify_destim() returned True")
            return True
        else:
            logger.debug("verify_destim() returned False")
            return False
        
    def tofile(self, fname, destim=None, hdrname=None, clobber=False):
        if not destim or not hdrname:
            self.hverify()
        self.writeto(fname, clobber=clobber)