#from .. mappings import DEGTORAD, RADTODEG
from hstwcs.mappings import DEGTORAD, RADTODEG
import numpy
from numpy import math 
from math import sin, sqrt, pow, cos, asin, atan2,pi

class MakeWCS(object):
    """
    Purpose
    =======
    Recompute WCS keywords based on PA_V3 and distortion model.
    
    Algorithm
    =========
    - update reference chip wcs
    - update extension wcs
    - update extension header
    
    """
    
    def __init__(self, owcs, refwcs):
        self.wcsobj = owcs
        self.refwcs = refwcs
        self.updateWCS()
        self.DOMakeWCS = 'COMPLETE'
        
    def updateWCS(self):    
        """
        recomputes the basic WCS kw 
        """
        backup_wcs = self.wcsobj.restore()
        ltvoff, offshift = self.getOffsets(backup_wcs)
            
        self.uprefwcs(ltvoff, offshift)
        self.upextwcs(ltvoff, offshift)
        self.updatehdr()
        
    def upextwcs(self, loff, offsh):
        """
        updates an extension wcs
        """  
        ltvoffx, ltvoffy = loff[0], loff[1]
        offshiftx, offshifty = offsh[0], offsh[1]
        backup_wcs = self.wcsobj.restore()
        ltv1 = self.wcsobj.ltv1
        ltv2 = self.wcsobj.ltv2  
        
        
        
        if ltv1 != 0. or ltv2 != 0.:
            offsetx = backup_wcs['CRPIX1'] - ltv1 - self.wcsobj.idcmodel.refpix['XREF']
            offsety = backup_wcs['CRPIX2'] - ltv2 - self.wcsobj.idcmodel.refpix['YREF']
            fx,fy = self.wcsobj.idcmodel.shift(self.wcsobj.idcmodel.cx,self.wcsobj.idcmodel.cy,offsetx,offsety)
        else:
            fx, fy = self.wcsobj.idcmodel.cx, self.wcsobj.idcmodel.cy
        
        ridcmodel = self.refwcs.idcmodel
        v2 = self.wcsobj.idcmodel.refpix['V2REF']
        v3 = self.wcsobj.idcmodel.refpix['V3REF']
        v2ref = self.refwcs.idcmodel.refpix['V2REF']

        v3ref = self.refwcs.idcmodel.refpix['V3REF']
        R_scale = self.refwcs.idcmodel.refpix['PSCALE']/3600.0
        off = sqrt((v2-v2ref)**2 + (v3-v3ref)**2)/(R_scale*3600.0)

        # Here we must include the PARITY
        if v3 == v3ref:
           theta=0.0
        else:
           theta = atan2(self.wcsobj.parity[0][0]*(v2-v2ref),self.wcsobj.parity[1][1]*(v3-v3ref))
        
        if self.refwcs.idcmodel.refpix['THETA']: theta += self.refwcs.idcmodel.refpix['THETA']*pi/180.0

        dX=(off*sin(theta)) + offshiftx
        dY=(off*cos(theta)) + offshifty
        px = numpy.array([[dX,dY]])
        newcrval = self.refwcs.pixel2world_fits(px)[0]
        newcrpix = numpy.array([self.wcsobj.idcmodel.refpix['XREF'] + ltvoffx, 
                                self.wcsobj.idcmodel.refpix['YREF'] + ltvoffy])
        
        self.wcsobj.crval = newcrval
        self.wcsobj.crpix = newcrpix
        
        # Account for subarray offset
        # Angle of chip relative to chip
        if self.wcsobj.idcmodel.refpix['THETA']:
           dtheta = self.wcsobj.idcmodel.refpix['THETA'] - self.refwcs.idcmodel.refpix['THETA']
        else:
           dtheta = 0.0
        
        
        # Create a small vector, in reference image pixel scale
        # There is no parity effect here ???
        
        delXX=fx[1,1]/R_scale/3600.
        delYX=fy[1,1]/R_scale/3600.
        delXY=fx[1,0]/R_scale/3600.
        delYY=fy[1,0]/R_scale/3600.
        
        
        # Convert to radians
        rr=dtheta*pi/180.0
        
        # Rotate the vectors
        dXX = cos(rr)*delXX - sin(rr)*delYX
        dYX = sin(rr)*delXX + cos(rr)*delYX

        dXY = cos(rr)*delXY - sin(rr)*delYY
        dYY = sin(rr)*delXY + cos(rr)*delYY
        px = numpy.array([[dX+dXX,dY+dYX]])
        
        # Transform to sky coordinates
        
        wc = self.refwcs.pixel2world_fits(px)
        
        a = wc[0,0]
        b = wc[0,1]
        px = numpy.array([[dX+dXY,dY+dYY]])
        
        wc = self.refwcs.pixel2world_fits(px)
        c = wc[0,0]
        d = wc[0,1]
        
        # Calculate the new CDs and convert to degrees
        cd11 = diff_angles(a,newcrval[0])*cos(newcrval[1]*pi/180.0)
        cd12 = diff_angles(c,newcrval[0])*cos(newcrval[1]*pi/180.0)
        cd21 = diff_angles(b,newcrval[1])
        cd22 = diff_angles(d,newcrval[1])
        
        cd = numpy.array([[cd11, cd12], [cd21, cd22]])
        self.wcsobj.cd = cd
        
    def uprefwcs(self, loff, offsh):
        """
        Updates the reference chip
        """
        ltvoffx, ltvoffy = loff[0], loff[1]
        offshift = offsh
        backup_refwcs = self.refwcs.restore()
        dec = backup_refwcs['CRVAL2']
    
        # Get an approximate reference position on the sky
        rref = numpy.array([[self.refwcs.idcmodel.refpix['XREF']+ltvoffx, 
                            self.refwcs.idcmodel.refpix['YREF']+ltvoffy]])
        
        
        crval = self.refwcs.pixel2world_fits(rref)
        
        # Convert the PA_V3 orientation to the orientation at the aperture
        # This is for the reference chip only - we use this for the
        # reference tangent plane definition
        # It has the same orientation as the reference chip
        pv = troll(self.wcsobj.pav3,dec,self.refwcs.idcmodel.refpix['V2REF'],
                    self.refwcs.idcmodel.refpix['V3REF'])
        # Add the chip rotation angle
        if self.refwcs.idcmodel.refpix['THETA']:
            pv += self.refwcs.idcmodel.refpix['THETA']
            
        
        # Set values for the rest of the reference WCS
        self.refwcs.crval = crval[0]
        self.refwcs.crpix = numpy.array([0.0,0.0])+offsh
        parity = self.refwcs.parity
        R_scale = self.refwcs.idcmodel.refpix['PSCALE']/3600.0
        cd11 = parity[0][0] *  cos(pv*pi/180.0)*R_scale
        cd12 = parity[0][0] * -sin(pv*pi/180.0)*R_scale
        cd21 = parity[1][1] *  sin(pv*pi/180.0)*R_scale
        cd22 = parity[1][1] *  cos(pv*pi/180.0)*R_scale
        rcd = numpy.array([[cd11, cd12], [cd21, cd22]])
        self.refwcs.cd = rcd
        self.refwcs.set()
    
    
        
    def updatehdr(self):
        """
        Keywords to update:
        CD1_1, CD1_2, CD2_1, CD2_2, CRPIX1, CRPIX2, CRVAL1, CRVAL2
        """
        
        self.wcsobj.hdr.update('CD1_1', self.wcsobj.cd[0,0])
        self.wcsobj.hdr.update('CD1_2', self.wcsobj.cd[0,1])
        self.wcsobj.hdr.update('CD2_1', self.wcsobj.cd[1,0])
        self.wcsobj.hdr.update('CD2_2', self.wcsobj.cd[1,1])
        self.wcsobj.hdr.update('CRVAL1', self.wcsobj.crval[0])
        self.wcsobj.hdr.update('CRVAL2', self.wcsobj.crval[1])
        self.wcsobj.hdr.update('CRPIX1', self.wcsobj.crpix[0])
        self.wcsobj.hdr.update('CRPIX2', self.wcsobj.crpix[1])
        self.wcsobj.hdr.update('ORIENTAT', self.wcsobj.orientat)
        
    def getOffsets(self, backup_wcs):
        ltv1 = self.wcsobj.ltv1
        ltv2 = self.wcsobj.ltv2
        
        offsetx = backup_wcs['CRPIX1'] - ltv1 - self.wcsobj.idcmodel.refpix['XREF']
        offsety = backup_wcs['CRPIX2'] - ltv2 - self.wcsobj.idcmodel.refpix['YREF']
        
        shiftx = self.wcsobj.idcmodel.refpix['XREF'] + ltv1
        shifty = self.wcsobj.idcmodel.refpix['YREF'] + ltv2
        if ltv1 != 0. or ltv2 != 0.:
            ltvoffx = ltv1 + offsetx
            ltvoffy = ltv2 + offsety
            offshiftx = offsetx + shiftx
            offshifty = offsety + shifty
        else:
            ltvoffx = 0.
            ltvoffy = 0.
            offshiftx = 0.
            offshifty = 0.
            
        ltvoff = numpy.array([ltvoffx, ltvoffy])
        offshift = numpy.array([offshiftx, offshifty])
        return ltvoff, offshift
    
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

def troll(roll, dec, v2, v3):
    """ Computes the roll angle at the target position based on:
            the roll angle at the V1 axis(roll),
            the dec of the target(dec), and
            the V2/V3 position of the aperture (v2,v3) in arcseconds.

        Based on the algorithm provided by Colin Cox that is used in
        Generic Conversion at STScI.
    """
    # Convert all angles to radians
    _roll = DEGTORAD(roll)
    _dec = DEGTORAD(dec)
    _v2 = DEGTORAD(v2 / 3600.)
    _v3 = DEGTORAD(v3 / 3600.)

    # compute components
    sin_rho = sqrt((pow(sin(_v2),2)+pow(sin(_v3),2)) - (pow(sin(_v2),2)*pow(sin(_v3),2)))
    rho = asin(sin_rho)
    beta = asin(sin(_v3)/sin_rho)
    if _v2 < 0: beta = pi - beta
    gamma = asin(sin(_v2)/sin_rho)
    if _v3 < 0: gamma = pi - gamma
    A = pi/2. + _roll - beta
    B = atan2( sin(A)*cos(_dec), (sin(_dec)*sin_rho - cos(_dec)*cos(rho)*cos(A)))

    # compute final value
    troll = RADTODEG(pi - (gamma+B))

    return troll

