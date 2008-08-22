#from .. mappings import DEGTORAD, RADTODEG
from hstwcs.mappings import DEGTORAD, RADTODEG
import numpy
from numpy import math 
from math import sin, sqrt, pow, cos, asin, atan2,pi
from hstwcs import utils

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
 
    def updateWCS(cls, wcs, refwcs):    
        """
        recomputes the basic WCS kw 
        """
        ltvoff, offshift = cls.getOffsets(wcs)
            
        cls.uprefwcs(wcs, refwcs, ltvoff, offshift)
        cls.upextwcs(wcs, refwcs, ltvoff, offshift)
        
        kw2update = {'CD1_1': wcs.cd[0,0],
                    'CD1_2': wcs.cd[0,1],
                    'CD2_1': wcs.cd[1,0],
                    'CD2_2': wcs.cd[1,1],
                    'CRVAL1': wcs.crval[0],
                    'CRVAL2': wcs.crval[1],
                    'CRPIX1': wcs.crpix[0],
                    'CRPIX2': wcs.crpix[1],
            }
        return kw2update
    
    updateWCS = classmethod(updateWCS)
    
    def upextwcs(self, wcs, refwcs, loff, offsh):
        """
        updates an extension wcs
        """  
        ltvoffx, ltvoffy = loff[0], loff[1]
        offshiftx, offshifty = offsh[0], offsh[1]
        ltv1 = wcs.ltv1
        ltv2 = wcs.ltv2
        
        if ltv1 != 0. or ltv2 != 0.:
            offsetx = backup_wcs['CRPIX1'] - ltv1 - wcs.idcmodel.refpix['XREF']
            offsety = backup_wcs['CRPIX2'] - ltv2 - wcs.idcmodel.refpix['YREF']
            fx,fy = wcs.idcmodel.shift(wcs.idcmodel.cx,wcs.idcmodel.cy,offsetx,offsety)
        else:
            fx, fy = wcs.idcmodel.cx, wcs.idcmodel.cy
        
        ridcmodel = refwcs.idcmodel
        v2 = wcs.idcmodel.refpix['V2REF']
        v3 = wcs.idcmodel.refpix['V3REF']
        v2ref = refwcs.idcmodel.refpix['V2REF']

        v3ref = refwcs.idcmodel.refpix['V3REF']
        R_scale = refwcs.idcmodel.refpix['PSCALE']/3600.0
        off = sqrt((v2-v2ref)**2 + (v3-v3ref)**2)/(R_scale*3600.0)
        if v3 == v3ref:
           theta=0.0
        else:
           theta = atan2(wcs.parity[0][0]*(v2-v2ref), wcs.parity[1][1]*(v3-v3ref))
        
        if refwcs.idcmodel.refpix['THETA']: theta += refwcs.idcmodel.refpix['THETA']*pi/180.0

        dX=(off*sin(theta)) + offshiftx
        dY=(off*cos(theta)) + offshifty
        px = numpy.array([[dX,dY]])
        newcrval = refwcs.pixel2world_fits(px)[0]
        newcrpix = numpy.array([wcs.idcmodel.refpix['XREF'] + ltvoffx, 
                                wcs.idcmodel.refpix['YREF'] + ltvoffy])
        wcs.crval = newcrval
        wcs.crpix = newcrpix
        
        # Account for subarray offset
        # Angle of chip relative to chip
        if wcs.idcmodel.refpix['THETA']:
           dtheta = wcs.idcmodel.refpix['THETA'] - refwcs.idcmodel.refpix['THETA']
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
        wc = refwcs.pixel2world_fits(px)
        
        a = wc[0,0]
        b = wc[0,1]
        px = numpy.array([[dX+dXY,dY+dYY]])
        
        wc = refwcs.pixel2world_fits(px)
        c = wc[0,0]
        d = wc[0,1]

        # Calculate the new CDs and convert to degrees
        cd11 = utils.diff_angles(a,newcrval[0])*cos(newcrval[1]*pi/180.0)
        cd12 = utils.diff_angles(c,newcrval[0])*cos(newcrval[1]*pi/180.0)
        cd21 = utils.diff_angles(b,newcrval[1])
        cd22 = utils.diff_angles(d,newcrval[1])
        cd = numpy.array([[cd11, cd12], [cd21, cd22]])
        wcs.cd = cd
        
    upextwcs = classmethod(upextwcs)
        
    def uprefwcs(self, wcs, refwcs, loff, offsh):
        """
        Updates the reference chip
        """
        ltvoffx, ltvoffy = loff[0], loff[1]
        offshift = offsh
        dec = refwcs.crval[1]
        # Get an approximate reference position on the sky
        rref = numpy.array([[refwcs.idcmodel.refpix['XREF']+ltvoffx, 
                            refwcs.idcmodel.refpix['YREF']+ltvoffy]])
        
        crval = refwcs.pixel2world_fits(rref)
        # Convert the PA_V3 orientation to the orientation at the aperture
        # This is for the reference chip only - we use this for the
        # reference tangent plane definition
        # It has the same orientation as the reference chip
        pv = troll(wcs.pav3,dec,refwcs.idcmodel.refpix['V2REF'],
                    refwcs.idcmodel.refpix['V3REF'])
        # Add the chip rotation angle
        if refwcs.idcmodel.refpix['THETA']:
            pv += refwcs.idcmodel.refpix['THETA']
            
        
        # Set values for the rest of the reference WCS
        refwcs.crval = crval[0]
        refwcs.crpix = numpy.array([0.0,0.0])+offsh
        parity = refwcs.parity
        R_scale = refwcs.idcmodel.refpix['PSCALE']/3600.0
        cd11 = parity[0][0] *  cos(pv*pi/180.0)*R_scale
        cd12 = parity[0][0] * -sin(pv*pi/180.0)*R_scale
        cd21 = parity[1][1] *  sin(pv*pi/180.0)*R_scale
        cd22 = parity[1][1] *  cos(pv*pi/180.0)*R_scale
        rcd = numpy.array([[cd11, cd12], [cd21, cd22]])
        refwcs.cd = rcd
        refwcs.set()
        
    uprefwcs = classmethod(uprefwcs)
        

    def getOffsets(cls, wcs):
        ltv1 = wcs.ltv1
        ltv2 = wcs.ltv2
        
        offsetx = wcs.crpix[0] - ltv1 - wcs.idcmodel.refpix['XREF']
        offsety = wcs.crpix[1] - ltv2 - wcs.idcmodel.refpix['YREF']
        
        shiftx = wcs.idcmodel.refpix['XREF'] + ltv1
        shifty = wcs.idcmodel.refpix['YREF'] + ltv2
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
    
    getOffsets = classmethod(getOffsets)
    

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

