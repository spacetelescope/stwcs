#from .. pywcs.sip import SIP
from pywcs.sip import SIP
import pyfits
import instruments
#from .. distortion import models
from hstwcs.distortion import models
import numpy as N

#from .. mappings import inst_mappings, ins_spec_kw, DEGTORAD, RADTODEG, basic_wcs
from hstwcs.mappings import inst_mappings, ins_spec_kw, DEGTORAD, RADTODEG
from hstwcs.mappings import basic_wcs, prim_hdr_kw

__docformat__ = 'restructuredtext'

class HSTWCS(SIP):
    """
    Purpose
    =======
    Create a WCS object based on the instrument.
    It has all basic WCS kw as attribbutes (set by pywcs).
    It also uses the primary and extension header to define 
    instrument specific attributes needed by the correction classes.
    """
    
    def __init__(self, hdr0, ehdr, fobj=None):
        """
        :Parameters:
        `hdr0`: Pyfits Header
                primary header
        `ehdr`: Pyfits Header
               extension header
        `fobj`: PyFITS HDUList object or None
                pyfits file object
        """
        SIP.__init__(self, ehdr, fobj=fobj)
        self.setHDR0kw(hdr0)
        self.detector = self.setDetector(hdr0)
        self.inst_kw = ins_spec_kw
        self.setInstrSpecKw(hdr0, ehdr)
        self.pscale = self.setPscale()
        self.orientat = self.setOrient()
        
        
       
    def setHDR0kw(self, primhdr):
        # Set attributes from kw defined in the primary header.
        self.instrument = primhdr.get('INSTRUME', None)
        self.offtab = primhdr.get('OFFTAB', None) 
        self.idctab = primhdr.get('IDCTAB', None)
        self.date_obs = primhdr.get('DATE-OBS', None)
        self.pav3 = primhdr.get('PA_V3', None)
        self.ra_targ = primhdr.get('RA_TARG', None)
        self.dec_targ = primhdr.get('DEC_TARG', None)
    
    
    def setDetector(self, header):
        # Set detector attribute for instuments which have more than one detector
        if self.instrument in ['ACS', 'WFC3']:
            return header.get('DETECTOR', None)
        else:
            return None
        
    def setInstrSpecKw(self, prim_hdr, ext_hdr):
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
        cd11 = self.cd[0][0]
        cd21 = self.cd[1][0]
        return N.sqrt(N.power(cd11,2)+N.power(cd21,2)) * 3600.
    
    def setOrient(self):
        # Recompute ORIENTAT
        cd12 = self.cd[0][1]
        cd22 = self.cd[1][1]
        return RADTODEG(N.arctan2(cd12,cd22))
    
    def verifyPa_V3(self):
        """make sure PA_V3 is populated"""    
        
    def verifyKw(self):
        """verify that all required kw have meaningful values"""
            
    def readModel(self, header):
        """
        Purpose
        =======
        Read distortion model from idc table.
        Save some of the information as kw needed later by pydrizzle
        
        list of kw - to be revised later
        """
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
        ext_hdr.update('IDCXSIZE', self.idcmodel.refpix['XSIZE'])
        ext_hdr.update('IDCYSIZE', self.idcmodel.refpix['YSIZE'])
        ext_hdr.update('IDCXDELT', self.idcmodel.refpix['XDELTA'])
        ext_hdr.update('IDCYDELT', self.idcmodel.refpix['YDELTA'])
        ext_hdr.update('CENTERED', self.idcmodel.refpix['centered'])
        
    
    