#from .. pywcs.sip import SIP
from pywcs import WCS
import pyfits
import instruments
#from .. distortion import models
from hstwcs.distortion import models
import numpy as N
from pytools import fileutil

#from .. mappings import inst_mappings, ins_spec_kw, DEGTORAD, RADTODEG, basic_wcs
from hstwcs.mappings import inst_mappings, ins_spec_kw, DEGTORAD, RADTODEG
from hstwcs.mappings import basic_wcs, prim_hdr_kw

__docformat__ = 'restructuredtext'

class HSTWCS(WCS):
    """
    Purpose
    =======
    Create a WCS object based on the instrument.
    It has all basic WCS kw as attribbutes (set by pywcs).
    It also uses the primary and extension header to define 
    instrument specific attributes needed by the correction classes.
    """
    
    def __init__(self, hdr0=None, ehdr=None, fobj=None, fname=None, instrument=None):
        """
        :Parameters:
        `fname`: string
                file/extension name
                filename[EXTNAME,EXTNUM]
                filename[extension_number]
        `hdr0`: Pyfits Header
                primary header
        `ehdr`: Pyfits Header
               extension header
        `fobj`: PyFITS HDUList object or None
                pyfits file object
        `instrument`: string
                one of 'ACS', 'NICMOS', 'WFPC2', 'STIS', 'WFC3' 
        """
        
        self.inst_kw = ins_spec_kw
        
        if fname == None and hdr0 == None and ehdr == None and instrument != None:
            # create a default HSTWCS object based on instrument only
            self.instrument = instrument
            self.setInstrSpecKw()
        elif fname != None:
            # create an HSTWCS object from a filename
            filename, extname = fileutil.parseFilename(fname)
            self.filename = filename
            self.extname = extname
            if extname == None:
                #data may be in the primary array
                ehdr = pyfits.getheader(filename)
                hdr0 = None
            else:
                ext = fileutil.parseExtn(extname)
                hdr0 = pyfits.getheader(filename)
                try:
                    ehdr = pyfits.getheader(filename, ext=ext)
                except:
                    print 'Unable to get extension header based on filename %s. /n' % fname
        elif fname == None and  instrument == None:
            # hdr0 may be None, a WCS object will still be created
            if hdr0 == None and ehdr == None:
                print 'Not enough information to create a WCS object\n'
                self.help()
                return
            if hdr0 != None:
                assert isinstance (hdr0, pyfits.Header)
            if ehdr != None:
                assert isinstance (ehdr, pyfits.Header)
                
        WCS.__init__(self, ehdr, fobj=fobj)    
        self.setHDR0kw(hdr0, ehdr)
        #self.detector = self.setDetector(hdr0)
        
        self.setInstrSpecKw(hdr0, ehdr)
        self.pscale = self.setPscale()
        self.orientat = self.setOrient()
        self.readIDCCoeffs(ehdr)
        # if input was not a file name, try to get it from the primary header
        self.filename = hdr0.get('FILENAME', "")
        extname = ehdr.get('EXTNAME', "")
        extnum = ehdr.get('EXTVER', None)
        self.extname = (extname, extnum)
            
        
        
    def setHDR0kw(self, primhdr, ehdr):
        if primhdr == None:
            # we are given only an extension header
            header = ehdr
        elif ehdr == None:
            header = primhdr
        else:
            hcards = primhdr.ascardlist()
            hcards.extend(ehdr.ascardlist())
            header = pyfits.Header(cards = hcards)
        # Set attributes from kw defined in the primary header.
        self.instrument = header.get('INSTRUME', None)
        self.offtab = header.get('OFFTAB', None) 
        self.idctab = header.get('IDCTAB', None)
        self.date_obs = header.get('DATE-OBS', None)
        self.pav3 = header.get('PA_V3', None)
        self.ra_targ = header.get('RA_TARG', None)
        self.dec_targ = header.get('DEC_TARG', None)
        self.detector = header.get('DETECTOR', None)
        self.filename = header.get('FILENAME', "")
    """
    def setDetector(self, header):
        # Set detector attribute for instuments which have more than one detector
        if self.instrument in ['ACS', 'WFC3']:
            return header.get('DETECTOR', None)
        else:
            return None
    """
    def readIDCCoeffs(self, header):
        """
        Reads in first order IDCTAB coefficients if present in the header
        """
        coeffs = ['ocx10', 'ocx11', 'ocy10', 'ocy11', 'idcscale']
        for c in coeffs:
            self.__setattr__(c, header.get(c, None))
            
    def setInstrSpecKw(self, prim_hdr=None, ext_hdr=None):
        # Based on the instrument kw creates an instalnce of an instrument WCS class
        # and sets attributes from instrument specific kw 
        if self.instrument in inst_mappings.keys():
            inst_kl = inst_mappings[self.instrument]
            inst_kl = instruments.__dict__[inst_kl]
            insobj = inst_kl(prim_hdr, ext_hdr)
            for key in self.inst_kw:
                self.__setattr__(key, insobj.__getattribute__(key))
        else:
            raise KeyError, "Unsupported instrument - %s" %self.instrument
    
    def setPscale(self):
        # Calculates the plate scale from the cd matrix
        # this may not be needed now and shoufl probably be done after all 
        # corrections
        cd11 = self.wcs.cd[0][0]
        cd21 = self.wcs.cd[1][0]
        return N.sqrt(N.power(cd11,2)+N.power(cd21,2)) * 3600.
    
    def setOrient(self):
        # Recompute ORIENTAT
        cd12 = self.wcs.cd[0][1]
        cd22 = self.wcs.cd[1][1]
        return RADTODEG(N.arctan2(cd12,cd22))
    
    def updatePscale(self, pscale):
        """Given a plate scale, update the CD matrix"""
        old_pscale = self.pscale
        self.pscale = pscale
        self.wcs.cd = self.wcs.cd * pscale/old_pscale
        self.naxis1 = self.naxis1 * old_pscale/ pscale 
        self.naxis2 = self.naxis2 * old_pscale/ pscale 
        self.wcs.crpix = self.wcs.crpix *old_pscale/pscale
        
    def updateOrient(self, orient):
        """Given n angle update the CD matrix"""
        if self.orientat == orient:
            return
        old_orient = self.orientat
        self.orientat = orient
        angle = fileutil.DEGTORAD(orient)
        cd11 = -N.cos(angle)
        cd12 = N.sin(angle)
        cd21 = cd12
        cd22 = -cd11
        cdmat = N.array([[cd11, cd12],[cd21,cd22]])
        self.wcs.cd = cdmat * self.pscale/3600
        
    def verifyPa_V3(self):
        """make sure PA_V3 is populated"""    
        
    def verifyKw(self):
        """verify that all required kw have meaningful values"""
            
    def readModel(self, header):
        """
        Purpose
        =======
        Read distortion model from idc table.
        Save some of the information as kw needed for interpreting the distortion
        
        """
        if self.idctab == None or self.date_obs == None:
            print 'idctab or date_obs not available\n'
            self.idcmodel = None
            return
        if self.filter1 ==None and self.filter2 == None:
            'No filter information available\n'
            self.idcmodel = None
            return
        
        self.idcmodel = models.IDCModel(self.idctab,
                    chip=self.chip, direction='forward', date=self.date_obs,
                    filter1=self.filter1, filter2=self.filter2, 
                    offtab=self.offtab, binned=self.binned)
        
        
        self.updatehdr(header)
                   
    
    def updatehdr(self, ext_hdr, newkeywords=None):
        #kw2add : OCX10, OCX11, OCY10, OCY11 
        # record the model in the header for use by pydrizzle
        ext_hdr.update('OCX10', self.idcmodel.cx[1,0])
        ext_hdr.update('OCX11', self.idcmodel.cx[1,1])
        ext_hdr.update('OCY10', self.idcmodel.cy[1,0])
        ext_hdr.update('OCY11', self.idcmodel.cy[1,1])
        ext_hdr.update('IDCSCALE', self.idcmodel.refpix['PSCALE'])
        ext_hdr.update('IDCTHETA', self.idcmodel.refpix['THETA'])
        ext_hdr.update('IDCXREF', self.idcmodel.refpix['XREF'])
        ext_hdr.update('IDCYREF', self.idcmodel.refpix['YREF'])
        ext_hdr.update('IDCV2REF', self.idcmodel.refpix['V2REF'])
        ext_hdr.update('IDCV3REF', self.idcmodel.refpix['V3REF'])
        
    
    def help(self):
        print 'How to create an HSTWCS object:\n\n'
        print """ \
        1. Create an HSTWCS object using pyfits.Header objects. \n
        Example:\n
        hdr0 = pyfits.getheader('j9irw4b1q_flt.fits', ext=0)\n
        hdr1 = pyfits.getheader('j9irw4b1q_flt.fits', ext=1)\n
        f = pyfits.open('j9irw4b1q_flt.fits')\n
        f is required only if a lookup table distortion is available.\n
        w = wcsutil.HSTWCS(hdr0, hdr1,f)\n\n
        
        2. Create an HSTWCS object using a qualified file name. \n
        Example:\n
        w = wcsutil.HSTWCS('j9irw4b1q_flt.fits[sci,1]')\n\n
        """
        
    