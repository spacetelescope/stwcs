import pyfits
import numpy as N
#from .. mappings import ins_spec_kw
from mappings import ins_spec_kw, prim_hdr_kw
    
class InstrWCS(object):
    """
    A base class for instrument specific keyword definition.
    It prvides a default implementation (modeled by ACS) for
    all set_kw methods.
    """
    def __init__(self, hdr0=None, hdr=None):
        self.exthdr = hdr   
        self.primhdr = hdr0
        
    def set_ins_spec_kw(self):
        """
        This method MUST call all set_kw methods.
        There should be a set_kw method for all kw listed in 
        mappings.ins_spec_kw
        """
        self.set_detector()
        self.set_filter1()
        self.set_filter2()
        self.set_vafactor()
        self.set_naxis1()
        self.set_naxis2()
        self.set_ltv1()
        self.set_ltv2()
        self.set_binned()
        self.set_chip()
        self.set_parity()
        
        
    def set_filter1(self):
        try:
            self.filter1 = self.primhdr.get('FILTER1', None)
        except:
            self.filter1 = None
    def set_filter2(self):
        try:
            self.filter2 = self.primhdr.get('FILTER2', None) 
        except:
            self.filter2 = None
    def set_vafactor(self):
        try:
            self.vafactor = self.exthdr.get('vafactor', 1) 
        except:
            self.vafactor = 1.
    def set_naxis1(self):
        try:
            self.naxis1 = self.exthdr.get('naxis1', None)
        except:
            self.naxis1 = None
    def set_naxis2(self):
        try:
            self.naxis2 = self.exthdr.get('naxis2', None)
        except:
            self.naxis2 = None
    def set_ltv1(self):
        try:
            self.ltv1 = self.exthdr.get('ltv1', 0.0)
        except:
            self.ltv1 = 0.0
        
    def set_ltv2(self):
        try:
            self.ltv2 = self.exthdr.get('ltv2', 0.0)
        except:
            self.ltv2 = 0.0
        
    def set_binned(self):
        try:
            self.binned = self.exthdr.get('BINAXIS1', 1)
        except:
            self.binned = 1
    def set_chip(self):
        try:
            self.chip = self.exthdr.get('CCDCHIP', 1)
        except:
            self.chip = 1
            
    def set_parity(self):
        self.parity = [[1.0,0.0],[0.0,-1.0]]
    
    def set_detector(self):
        # each instrument has a different kw for detector and it can be 
        # in a different header, so this is to be handled by the instrument classes 
        pass
        
class ACSWCS(InstrWCS):
    """
    get instrument specific kw    
    """
    
    def __init__(self, hdr0, hdr):
        self.primhdr = hdr0
        self.exthdr = hdr
        InstrWCS.__init__(self,hdr0, hdr)
        self.set_ins_spec_kw()
        
    def set_detector(self):
        try:
            self.detector = self.primhdr['DETECTOR']  
        except KeyError:
            print 'ERROR: Detector kw not found.\n'
            raise
        except TypeError:
            #this is the case of creating a default HSTWCS object by 
            #providing 'instrument' and 'detector'
            pass
    
    def set_parity(self):
        parity = {'WFC':[[1.0,0.0],[0.0,-1.0]],
                'HRC':[[-1.0,0.0],[0.0,1.0]],
                'SBC':[[-1.0,0.0],[0.0,1.0]]}
        
        if self.detector not in parity.keys():
            parity = InstrWCS.set_parity(self)
        else:
            self.parity = parity[self.detector]
    
        
class WFPC2WCS(InstrWCS):   


    def __init__(self, hdr0, hdr):
        self.primhdr = hdr0
        self.exthdr = hdr
        InstrWCS.__init__(self,hdr0, hdr)
        self.set_ins_spec_kw()
    
    def set_filter1(self):
        self.filter1 = self.primhdr.get('FILTNAM1', None)
        if self.filter1 == " " or self.filter1 == None:
            self.filter1 = 'CLEAR1'

    def set_filter2(self):
        self.filter2 = self.primhdr.get('FILTNAM2', None)
        if self.filter2 == " " or self.filter2 == None:
            self.filter2 = 'CLEAR2'
            
    
    def set_binned(self):
        mode = self.primhdr.get('MODE', 1)
        if mode == 'FULL':
            self.binned = 1
        elif mode == 'AREA':
            self.binned = 2

    def set_chip(self):
        self.chip = self.exthdr.get('DETECTOR', 1)
    
    def set_parity(self):
        self.parity = [[-1.0,0.],[0.,1.0]]
        
    def set_detector(self):
        try:
            self.detector = self.exthdr['DETECTOR']
        except KeyError:
            print 'ERROR: Detector kw not found.\n'
            raise
    

class WFC3WCS(InstrWCS):
    """
    Create a WFC3 detector specific class
    """
    
    def __init__(self, hdr0, hdr):
        self.primhdr = hdr0
        self.exthdr = hdr
        InstrWCS.__init__(self,hdr0, hdr)
        self.set_ins_spec_kw()
    
    def set_detector(self):
        try:
            self.detector = self.primhdr['DETECTOR']  
        except KeyError:
            print 'ERROR: Detector kw not found.\n'
            raise
        
    def set_filter1(self):
        self.filter1 = self.primhdr.get('FILTER', None)
        if self.filter1 == " " or self.filter1 == None:
            self.filter1 = 'CLEAR'
    
    def set_filter2(self):
        #Nicmos idc tables do not allow 2 filters.
        self.filter2 = 'CLEAR'
        
    def set_parity(self):
        parity = {'UVIS':[[-1.0,0.0],[0.0,1.0]], 
          'IR':[[-1.0,0.0],[0.0,1.0]]}
        
        if self.detector not in parity.keys():
            parity = InstrWCS.set_parity(self)
        else:
            self.parity = parity[self.detector]
            
class NICMOSWCS(InstrWCS):
    """
    Create a NICMOS specific class
    """
    
    def __init__(self, hdr0, hdr):
        self.primhdr = hdr0
        self.exthdr = hdr
        InstrWCS.__init__(self,hdr0, hdr)
        self.set_ins_spec_kw()
    
    def set_parity(self):
        self.parity = [[-1.0,0.],[0.,1.0]]
        
    def set_filter1(self):
        self.filter1 = self.primhdr.get('FILTER', None)
        if self.filter1 == " " or self.filter1 == None:
            self.filter1 = 'CLEAR'

    def set_filter2(self):
        #Nicmos idc tables do not allow 2 filters.
        self.filter2 = 'CLEAR'
        """
        self.filter2 = self.primhdr.get('FILTER2', None)
        if self.filter2 == " " or self.filter2 == None:
            self.filter2 = 'CLEAR2'
        """
    def set_chip(self):
        self.chip = self.detector
            
    def set_detector(self):
        try:
            self.detector = self.primhdr['CAMERA']  
        except KeyError:
            print 'ERROR: Detector kw not found.\n'
            raise
    
class STISWCS(InstrWCS):
    """
    Create a NICMOS specific class
    """
    
    def __init__(self, hdr0, hdr):
        self.primhdr = hdr0
        self.exthdr = hdr
        InstrWCS.__init__(self,hdr0, hdr)
        self.set_ins_spec_kw()
    
    def set_parity(self):
        self.parity = [[-1.0,0.],[0.,1.0]]
        
    def set_filter1(self):
        self.filter1 = self.exthdr.get('OPT_ELEM', None)
        if self.filter1 == " " or self.filter1 == None:
            self.filter1 = 'CLEAR1'

    def set_filter2(self):
        self.filter2 = self.exthdr.get('FILTER', None)
        if self.filter2 == " " or self.filter2 == None:
            self.filter2 = 'CLEAR2'
            
    def set_detector(self):
        try:
            self.detector = self.primhdr['DETECTOR']  
        except KeyError:
            print 'ERROR: Detector kw not found.\n'
            raise