import datetime
import numpy
from numpy import linalg
from pytools import fileutil
from hstwcs.utils import diff_angles
import makewcs, dgeo

MakeWCS = makewcs.MakeWCS
DGEO = dgeo.DGEO

class TDDCorr(object):
    """
    Purpose
    =======
    Apply time dependent distortion correction to SIP coefficients and basic
    WCS keywords. It is applicable only to ACS/WFC data.
    
    Ref: Jay Anderson, ACS ISR 2007-08, Variation of the Distortion 
    Solution of the WFC 
    
    :Parameters:
    `owcs`: HSTWCS object
            An extension HSTWCS object to be modified
    `refwcs`: HSTWCS object
             A reference HSTWCS object
    """
    
    def updateWCS(cls, wcs, refwcs):
        """
        - Calculates alpha and beta for ACS/WFC data.
        - Writes 2 new kw to the extension header: TDDALPHA and TDDBETA
        - sets TDDCORR to COMPLETE
        """
        
        alpha, beta = cls.set_alpha_beta(wcs)
        cls.apply_tdd2idc(wcs, alpha, beta)
        newkw = {'TDDALPHA': alpha, 'TDDBETA':beta, 'TDDCORR': 'COMPLETE'}
        
        return newkw
    updateWCS = classmethod(updateWCS)
    
    def apply_tdd2idc(cls,wcs, alpha, beta):
        """
        Applies TDD to the idctab coefficients of a ACS/WFC observation.
        This should be always the first correction.
        """
        theta_v2v3 = 2.234529
        idctheta = theta_v2v3
        idcrad = fileutil.DEGTORAD(idctheta)
        mrotp = fileutil.buildRotMatrix(idctheta)
        mrotn = fileutil.buildRotMatrix(-idctheta)
        abmat = numpy.array([[beta,alpha],[alpha,beta]])
        tdd_mat = numpy.array([[1+(beta/2048.), alpha/2048.],[alpha/2048.,1-(beta/2048.)]],numpy.float64)
        abmat1 = numpy.dot(tdd_mat, mrotn)
        abmat2 = numpy.dot(mrotp,abmat1)
        xshape, yshape = wcs.idcmodel.cx.shape, wcs.idcmodel.cy.shape
        icxy = numpy.dot(abmat2,[wcs.idcmodel.cx.ravel(),wcs.idcmodel.cy.ravel()])
        wcs.idcmodel.cx = icxy[0]
        wcs.idcmodel.cy = icxy[1]
        wcs.idcmodel.cx.shape = xshape
        wcs.idcmodel.cy.shape = yshape
        wcs.set()
        
    apply_tdd2idc = classmethod(apply_tdd2idc)
        
    def set_alpha_beta(cls, wcs):
        """
        Compute the time dependent distortion skew terms
        default date of 2004.5 = 2004-7-1
        """
        dday = 2004.5
        year,month,day = wcs.date_obs.split('-')
        rdate = datetime.datetime(int(year),int(month),int(day))
        rday = float(rdate.strftime("%j"))/365.25 + rdate.year
        alpha = 0.095 + 0.090*(rday-dday)/2.5
        beta = -0.029 - 0.030*(rday-dday)/2.5
        
        return alpha, beta

    set_alpha_beta = classmethod(set_alpha_beta)
    
        
class VACorr(object):
    """
    Purpose
    =======
    Apply velocity aberation correction to WCS keywords.
    Modifies the CD matrix and CRVAL1/2

    """
    
    def updateWCS(cls, wcs, refwcs):
        if wcs.vafactor != 1:
            wcs.cd = wcs.cd * wcs.vafactor
            crval1 = wcs.crval[0]
            crval2 = wcs.crval[1]
            wcs.crval[0] = wcs.crval[0] + wcs.vafactor*diff_angles(wcs.crval[0], crval1) 
            wcs.crval[1] = wcs.crval[1] + wcs.vafactor*diff_angles(wcs.crval[1], crval2) 
            
            kw2update={'CD1_1': wcs.cd[0,0], 'CD1_2':wcs.cd[0,1], 
                    'CD2_1':wcs.cd[1,0], 'CD2_2':wcs.cd[1,1], 
                    'CRVAL1':wcs.crval[0], 'CRVAL2':wcs.crval[1]}
            wcs.set()
        else:
            pass
        return kw2update
    
    updateWCS = classmethod(updateWCS)

        
class CompSIP(object):
    """
    Purpose
    =======
    Compute SIP coefficients from idc table coefficients.
    """
    
    def updateWCS(cls, wcs, refwcs):
        kw2update = {}
        order = wcs.idcmodel.norder
        kw2update['A_ORDER'] = order
        kw2update['B_ORDER'] = order
        pscale = wcs.idcmodel.refpix['PSCALE']
        
        cx = wcs.idcmodel.cx
        cy = wcs.idcmodel.cy
        
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
        return kw2update
    
    updateWCS = classmethod(updateWCS)


