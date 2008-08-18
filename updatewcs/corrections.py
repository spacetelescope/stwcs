import datetime
import numpy
from numpy import linalg
from pytools import fileutil

from makewcs import MakeWCS
from dgeo import DGEO

class TDDCorr(object):
    """
    Purpose
    =======
    Apply time dependent distortion correction to SIP coefficients and basic
    WCS keywords. Atthi stime it is applicable only to ACS/WFC data.
    
    Ref: Jay Anderson, ACS ISR 2007-08, Variation of the Distortion 
    Solution of the WFC 
    
    """
    def __init__(self, owcs=None, refwcs=None):
        """
        :Parameters:
        `owcs`: HSTWCS object
                An extension HSTWCS object to be modified
        `refwcs`: HSTWCS object
                 A reference HSTWCS object
        """
        self.wcsobj = owcs
        if self.wcsobj.DOTDDCorr == 'PERFORM':
            self.updateWCS()
            self.wcsobj.DOTDDCorr = 'COMPLETE'
        else:
            pass

    def updateWCS(self):
        """
        - Calculates alpha and beta for ACS/WFC data.
        - Writes 2 new kw to the extension header: TDDALPHA and TDDBETA
        """
        if self.wcsobj.instrument == 'ACS' and self.wcsobj.detector == 'WFC':
            newkw = ['TDDALPHA', 'TDDBETA']
            self.set_alpha_beta()
            self.applyTDD()
            self.updatehdr(newkeywords=newkw)
            
        else:
            pass
    
    def applyTDD(self):
        """
        Applies TDD to the idctab coefficients of a ACS/WFC observation.
        This should be always the first correction.
        """
        theta_v2v3 = 2.234529
        idctheta = theta_v2v3
        idcrad = fileutil.DEGTORAD(idctheta)
        mrotp = fileutil.buildRotMatrix(idctheta)
        mrotn = fileutil.buildRotMatrix(-idctheta)
        abmat = numpy.array([[self.TDDBETA,self.TDDALPHA],[self.TDDALPHA,self.TDDBETA]])
        tdd_mat = numpy.array([[1+(self.TDDBETA/2048.), self.TDDALPHA/2048.],[self.TDDALPHA/2048.,1-(self.TDDBETA/2048.)]],numpy.float64)
        abmat1 = numpy.dot(tdd_mat, mrotn)
        abmat2 = numpy.dot(mrotp,abmat1)
        xshape, yshape = self.wcsobj.idcmodel.cx.shape, self.wcsobj.idcmodel.cy.shape
        icxy = numpy.dot(abmat2,[self.wcsobj.idcmodel.cx.ravel(),self.wcsobj.idcmodel.cy.ravel()])
        self.wcsobj.idcmodel.cx = icxy[0]
        self.wcsobj.idcmodel.cy = icxy[1]
        self.wcsobj.idcmodel.cx.shape = xshape
        self.wcsobj.idcmodel.cy.shape = yshape
        
    def set_alpha_beta(self):
        """
        Compute the time dependent distortion skew terms
        default date of 2004.5 = 2004-7-1
        """
        dday = 2004.5
        year,month,day = self.wcsobj.date_obs.split('-')
        rdate = datetime.datetime(int(year),int(month),int(day))
        rday = float(rdate.strftime("%j"))/365.25 + rdate.year
        alpha = 0.095 + 0.090*(rday-dday)/2.5
        beta = -0.029 - 0.030*(rday-dday)/2.5
        self.TDDALPHA = alpha
        self.TDDBETA = beta

    def updatehdr(self, newkeywords=None):
        
        for kw in newkeywords:
            self.wcsobj.hdr.update(kw, self.__getattribute__(kw))
        
        
class VACorr(object):
    """
    Purpose
    =======
    Apply velocity aberation correction to WCS keywords.
    Modifies the CD matrix and CRVAL1/2
    """
    def __init__(self, owcs=None, refwcs=None):
        self.wcsobj = owcs
        self.refwcs = refwcs
        if self.wcsobj.DOVACorr == 'PERFORM':
            self.updateWCS()
            self.wcsobj.DOVACorr == 'COMPLETE'
        else:
            pass
        
    def updateWCS(self):
        if self.wcsobj.vafactor != 1:
            self.wcsobj.cd = self.wcsobj.cd * self.wcsobj.vafactor
            
            #self.wcsobj.crval[0] = self.wcsobj.ra_targ + self.wcsobj.vafactor*(self.wcsobj.crval[0] - self.wcsobj.ra_targ) 
            #self.wcsobj.crval[1] = self.wcsobj.dec_targ + self.wcsobj.vafactor*(self.wcsobj.crval[1] - self.wcsobj.dec_targ) 
            ref_backup = self.refwcs.restore()
            crval1 = ref_backup['CRVAL1']
            crval2 = ref_backup['CRVAL2']
            self.wcsobj.crval[0] = crval1 + self.wcsobj.vafactor*diff_angles(self.wcsobj.crval[0], crval1) 
            self.wcsobj.crval[1] = crval2 + self.wcsobj.vafactor*diff_angles(self.wcsobj.crval[1], crval2) 

            kw2update={'CD1_1': self.wcsobj.cd[0,0], 'CD1_2':self.wcsobj.cd[0,1], 
                    'CD2_1':self.wcsobj.cd[1,0], 'CD2_2':self.wcsobj.cd[1,1], 
                    'CRVAL1':self.wcsobj.crval[0], 'CRVAL2':self.wcsobj.crval[1]}
            self.updatehdr(newkeywords=kw2update)
        else:
            pass
            
    def updatehdr(self, newkeywords=None):
        for kw in newkeywords:
            self.wcsobj.hdr.update(kw, newkeywords[kw])
        
class CompSIP(object):
    """
    Purpose
    =======
    Compute SIP coefficients from idc table coefficients.
        
    """
    def __init__(self, owcs, refwcs):
        self.wcsobj = owcs
        self.updateWCS()
        self.DOCOMPSIP = 'COMPLETE'
        
    def updateWCS(self):
        kw2update = {}
        order = self.wcsobj.idcmodel.norder
        kw2update['A_ORDER'] = order
        kw2update['B_ORDER'] = order
        pscale = self.wcsobj.idcmodel.refpix['PSCALE']
        
        cx = self.wcsobj.idcmodel.cx
        cy = self.wcsobj.idcmodel.cy
        
        matr = numpy.array([[cx[1,1],cx[1,0]], [cy[1,1],cy[1,0]]], dtype=numpy.float)
        imatr = linalg.inv(matr)
        akeys1 = numpy.zeros((order+1,order+1), dtype=numpy.float)
        bkeys1 = numpy.zeros((order+1,order+1), dtype=numpy.float)
        for n in range(order+1):
            for m in range(order+1):
                if n >= m and n>=2:
                    idcval = numpy.array([[cx[n,m]],[cy[n,m]]])
                    sipval = numpy.dot(imatr, idcval)
                    akeys1[m,n-m] = sipval[0]
                    bkeys1[m,n-m] = sipval[1]
                    Akey="A_%d_%d" % (m,n-m)
                    Bkey="B_%d_%d" % (m,n-m)
                    kw2update[Akey] = sipval[0,0]
                    kw2update[Bkey] = sipval[1,0]
        self.updatehdr(newkw=kw2update)
                    
    def updatehdr(self, newkw=None):
        if not newkw: return
        for kw in newkw.keys():
            self.wcsobj.hdr.update(kw, newkw[kw])
            
            
def diff_angles(a,b):
    """ 
    Perform angle subtraction a-b taking into account
    small-angle differences across 360degree line. 
    """
    
    diff = a - b

    if diff > 180.0:
       diff -= 360.0

    if diff < -180.0:
       diff += 360.0
    
    return diff
