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
    Create a wcs object based on the instrument.
    It has all basic WCS kw as attribbutes (set by pywcs).
    It also uses the primary and extension header to define 
    instrument specific attributes needed by the correction classes.
    """
    
    def __init__(self, hdr0, hdr):
        """
        :Parameters:
        `hdr0`: Pyfits Header
                primary header
        `hdr`: Pyfits Header
               extension header
        """
        self.hdr = hdr
        self.hdr0 = hdr0
        SIP.__init__(self, self.hdr)
        self.setHDR0kw(hdr0)
        self.detector = self.setDetector()
        self.inst_kw = ins_spec_kw
        self.setInstrSpecKw()
        self.pscale = self.setPscale()
        self.orientat = self.setOrient()
        
        
       
    def setHDR0kw(self, primhdr):
        #sets attributes from kw defined in the primary header
        self.instrument = primhdr.get('INSTRUME', None)
        self.offtab = primhdr.get('OFFTAB', None) 
        self.idctab = primhdr.get('IDCTAB', None)
        self.date_obs = primhdr.get('DATE-OBS', None)
        self.pav3 = primhdr.get('PA_V3', None)
        self.ra_targ = primhdr.get('RA_TARG', None)
        self.dec_targ = primhdr.get('DEC_TARG', None)
    
    
    def setDetector(self):
        #sets detector attribute for instuments which have more than one detector
        if self.instrument in ['ACS', 'WFC3']:
            return self.hdr0.get('DETECTOR', None)
    
    
    def setInstrSpecKw(self):
        #Based on the instrument kw creates an instalnce of an instrument WCS class
        #and sets attributes from instrument specific kw 
        if self.instrument in inst_mappings.keys():
            inst_kl = inst_mappings[self.instrument]
            inst_kl = instruments.__dict__[inst_kl]
            insobj = inst_kl(self.hdr0, self.hdr)
            for key in self.inst_kw:
                self.__setattr__(key, insobj.__getattribute__(key))
        else:
            raise KeyError, "Unsupported instrument - %s" %self.instrument
    
    def setPscale(self):
        #Calculates the plate scale from the cd matrix
        #this may not be needed now and shoufl probably be done after all 
        # corrections
        cd11 = self.cd[0][0]
        cd21 = self.cd[1][0]
        return N.sqrt(N.power(cd11,2)+N.power(cd21,2)) * 3600.
    
    def setOrient(self):
        #same as setPscale
        cd12 = self.cd[0][1]
        cd22 = self.cd[1][1]
        return RADTODEG(N.arctan2(cd12,cd22))
    
    def verifyPa_V3(self):
        """make sure PA_V3 is populated"""    
        
    def verifyKw(self):
        """verify that all required kw have meaningful values"""
            
    def readModel(self):
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
        
        
        self.updatehdr()
                   
    
    def updatehdr(self, newkeywords=None):
        #kw2add : OCX10, OCX11, OCY10, OCY11 
        # record the model in the header for use by pydrizzle
        self.hdr.update('OCX10', self.idcmodel.cx[1,0])
        self.hdr.update('OCX11', self.idcmodel.cx[1,1])
        self.hdr.update('OCY10', self.idcmodel.cy[1,0])
        self.hdr.update('OCY11', self.idcmodel.cy[1,1])
        self.hdr.update('IDCSCALE', self.idcmodel.refpix['PSCALE'])
        self.hdr.update('IDCTHETA', self.idcmodel.refpix['THETA'])
        self.hdr.update('IDCXREF', self.idcmodel.refpix['XREF'])
        self.hdr.update('IDCYREF', self.idcmodel.refpix['YREF'])
        self.hdr.update('IDCV2REF', self.idcmodel.refpix['V2REF'])
        self.hdr.update('IDCV3REF', self.idcmodel.refpix['V3REF'])
        self.hdr.update('IDCXSIZE', self.idcmodel.refpix['XSIZE'])
        self.hdr.update('IDCYSIZE', self.idcmodel.refpix['YSIZE'])
        self.hdr.update('IDCXDELT', self.idcmodel.refpix['XDELTA'])
        self.hdr.update('IDCYDELT', self.idcmodel.refpix['YDELTA'])
        self.hdr.update('CENTERED', self.idcmodel.refpix['centered'])
        
    
    def restore(self):
        """
        restores basic wcs keywords from archive
        this should be modified to always return a populated 
        dictionary, although the values may be None.
        """    

        backup = {}
        for k in basic_wcs:
            try:
                nkw = ('O'+k)[:7]
                #backup[k] = self.hdr.__getitem__('O'+k)
                backup[k] = self.hdr[nkw]
                #self.__setattr__(k, self.hdr.__getitem__('O'+k))
                self.__setattr__(k, self.hdr[nkw])
            except KeyError:
                print 'Keyword %s not found' % k
                
        return backup
    
    def archive_kw(self):
        """
        archive original WCS kw before recalculating them.
        """
        backup_kw = self.restore()
        if backup_kw != {}:
            return
        else:
            # keep the archived kw 8 characters long
            cmt = 'archived value'
            for kw in basic_wcs:
                nkw = 'O'+kw
                self.hdr.update(nkw[:7], self.hdr[kw], comment=cmt)
        

