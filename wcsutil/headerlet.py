from __future__ import division
import time
import pyfits
from hstwcs import HSTWCS
import altwcs
from pyfits import HDUList
from stwcs import __version__ as stwcsVersion
from mappings import basic_wcs


def createHeaderlet(fname, hdrname, destim=None, output=None):
    """
    creates a headerlet from a science observation
    """
    fmt="%Y-%m-%dT%H:%M:%S"
    phdukw = {'IDCTAB': True, 
            'NPOLFILE': True, 
            'D2IMFILE': True} 
    fobj = pyfits.open(fname)
    if destim is None:
        try:
            destim = fobj[0].header['ROOTNAME']
        except KeyError:
            raise ValueError, 'Please provide a value for the DESTIM keyword'
    if hdrname is None:
        raise ValueError, "Please provide a name for the headerlet, HDRNAME is a required parameter."
    if output is None:
        output = hdrname
    altkeys = altwcs.wcskeys(fobj[('SCI',1)].header)
    if 'O' in altkeys:
        altkeys.remove('O')
    numsci = countext(fname, 'SCI')
    hdul = pyfits.HDUList()
    phdu = pyfits.PrimaryHDU()
    phdu.header.update('DESTIM', destim, comment='Destination observation root name')
    phdu.header.update('HDRNAME',hdrname, 'Headerlet name')
    phdu.header.update('STWCSVER',stwcsVersion, comment='Version of STWCS that generated this file')
    phdu.header.update('DATE',time.strftime(fmt), comment='Date FITS file was generated')

    updateRefFiles(fobj[0].header.ascard, phdu.header.ascard)
    phdu.header.update(key='VAFACTOR', value=fobj[('SCI',1)].header.get('VAFACTOR', 1.))
    hdul.append(phdu)
    
    for e in range(1,numsci+1):
        hwcs = HSTWCS(fname,ext=('SCI',e))
        h = hwcs.wcs2header(sip2hdr=True).ascard
        for ak in altkeys:
            awcs = HSTWCS(fname,ext=('SCI',e), wcskey=ak)
            h.extend(awcs.wcs2header().ascard)
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
                print "Required keyword AXISCORR was not found in %s['SCI',%d]" %(fname,e)
                print "Unable to create headerlet ... quiting\n"
                raise
            
            try: h.append(fhdr['D2IMERR'])
            except KeyError: h.append(pyfits.Card(key='DPERROR', value=0, 
                                                comment='Maximum error of D2IMARR'))
            
        hdu = pyfits.ImageHDU(header=pyfits.Header(h))
        hdul.append(hdu)
    numwdvarr = countext(fname, 'WCSDVARR')
    numd2im = countext(fname, 'D2IMARR')
    for w in range(1, numwdvarr+1):
        hdul.append(fobj[('WCSDVARR',w)].copy())
    for d in range(1, numd2im+1):
        hdul.append(fobj[('D2IMARR',d)].copy())
    hdul.writeto(output,clobber=True)
    fobj.close()

def updateRefFiles(source, dest):
    #dest is destination extension header.ascard
    #source is source header.ascard
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
    
def mapFitsExt2HDUListInd(fname, extname):
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
    if isinstance(fname, str):
        f = pyfits.open(fname)
    else:
        f = fname
    numext = 0
    for hdu in f:
        if hdu.header.has_key('EXTNAME') and hdu.header['EXTNAME'] == extname:
            numext+=1
    return numext

def _cleanDestWCS(dest):
    #clean all WCSs in the destination file
    #dkeys = altwcs.wcskeys(pyfits.getheader(dest,ext=('SCI',1)))
    #Warning: Does not clean teh primary WCS
    fobj = pyfits.open(dest, mode='update')
    numext = len(fobj)
    
    for e in range(numext):
        removeD2IM(fobj[e])
        removeSIP(fobj[e])
        removeLUT(fobj[e])
        removePrimaryWCS(fobj,e) 
    removeAltWCS(fobj, ext=range(numext)) 
    numwdvarr = countext(dest,'WCSDVARR')
    numd2im = countext(dest,'D2IMARR')
    for i in range(1, numwdvarr+1):
        del fobj[('WCSDVARR',i)]
    for i in range(1, numd2im+1):
        del fobj[('D2IMARR',i)]
    fobj.close()
    
def removeSIP(ext):
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

def removeLUT(ext):
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

def removeD2IM(ext):
    d2imkeys = ['D2IMFILE', 'AXISCORR', 'D2IMEXT', 'D2IMERR']
    for k in d2imkeys:
        try:
            del ext.header[k]
        except KeyError:
            pass
        
def removeAltWCS(dest, ext):
    dkeys = altwcs.wcskeys(dest[('SCI',1)].header)
    for k in dkeys:
        altwcs.deleteWCS(dest, ext=ext, wcskey=k)

def removePrimaryWCS(dest, ext):
    naxis=dest[ext].header.ascard['NAXIS'].value
    for key in basic_wcs:
        for i in range(1,naxis+1):
            try:
                del dest[ext].header.ascard[key+str(i)]
            except KeyError: pass
    try:
        del dest[ext].header.ascard['WCSAXES']
    except KeyError: pass

class Headerlet(HDUList):
    def __init__(self, fname, wcskeys=[]):
        
        self.wcskeys = wcskeys
        if isinstance(fname,str):
            fobj = pyfits.open(fname)
            #fobj.close()
        elif isinstance(fname, list):
            fobj = fname
        else:
            raise ValueError, "Input must be a file name (string) or a pyfits file object (HDUList)."
        HDUList.__init__(self, fobj) 
        self.fname = self.filename()
        self.hdrname = self[0].header.get("HDRNAME","")
        self.stwcsver = stwcsVersion
        self.destim = self[0].header.get("DESTIM","")
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
            _cleanDestWCS(dest)
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
                for k in siphdr.keys():
                    if k not in ['XTENSION', 'BITPIX', 'NAXIS', 'PCOUNT', 
                                 'GCOUNT','EXTNAME', 'EXTVER', 'ORIGIN', 
                                 'INHERIT', 'DATE', 'IRAF-TLM']:
                        fhdr.insert(wind, siphdr[k])
                    else:
                        pass
            #! Always attach these extensions last. Otherwise thier headers may
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
            print "Observation %s cannot be updated with headerlet %s" %(dest, self.hdrname)
        
        
    def hverify(self):
        self.verify()
        assert(self[0].header.has_key('DESTIM') and self[0].header['DESTIM']!= " ")
        assert(self[0].header.has_key('HDRNAME') and self[0].header['HDRNAME']!= " ")
        assert(self[0].header.has_key('STWCSVER') and self[0].header['STWCSVER']!= " ")
        """
        npolfile, vafactor and d2imfile are optional.
        idctab may be optional too ...
        assert(self[0].header.has_key('IDCTAB') 
            and self[0].header['IDCTAB']!= " ")
        assert(self[0].header.has_key('NPOLFILE') 
            and self[0].header['NPOLFILE']!= " ")
        """
    
    def verify_dest(self, dest):
        """
        verifies that the headerlet can be applied to the observation
        requirements:
        - if self.destim = fname.rootname
        - if fname has the same file structure.
        - if idctab, npolfile, d2imfile are different for self and fname
        - if vafactor is different
        - 
        """
        dobj = pyfits.open(dest)
        try:
            droot = dobj[0].header['ROOTNAME']
        except KeyError:
            dobj.close()
            print "This headerlet cannot be applied to observation %s. " % destfile
            print "It can only be applied to observations with root name %s" % self.destim
            return False
        dnumsci = countext(dest, 'SCI')
        hnumsip = countext(self, 'SIPWCS')
        if dnumsci != hnumsip:
            print "This headerlet cannot be applied to observation %s. " % dest
            print "The observation and the headerlet have a different number of 'SCI' extensions."
            return False
        dobj.close()
        return True
        
    def tofile(self, fname, destim=None, hdrname=None, clobber=False):
        if not destim or not hdrname:
            self.hverify()
        self.writeto(fname, clobber=clobber)