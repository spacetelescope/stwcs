import numpy as np
import pywcs
import pyfits
from updatewcs import wcsutil
from numpy import sqrt, arctan2

def output_wcs(list_of_wcsobj, ref_wcs=None, outwcs=None):
    fra_dec = np.vstack([w.footprint for w in list_of_wcsobj])
    
    delta_fdec = (fra_dec[:,1].max()-fra_dec[:,1].min())
    crval2 = fra_dec[:,1].min()+ delta_fdec/2.
    delta_fra = (fra_dec[:,0].max()-fra_dec[:,0].min())
    min_ind = fra_dec[:,0].argmin()
    crval1 = (fra_dec[min_ind,0]+ (delta_fra/2.)*np.cos((fra_dec[min_ind,1]-crval2)*np.pi/180))
    
    crval = np.array([crval1,crval2])
    if outwcs is None:
        if ref_wcs == None:
            ref_wcs = list_of_wcsobj[0]
            
        outwcs = undistortWCS(ref_wcs)
        outwcs.wcs.crval = crval
        outwcs.pscale = sqrt(outwcs.wcs.cd[0,0]**2 + outwcs.wcs.cd[1,0]**2)*3600.
        outwcs.orientat = arctan2(outwcs.wcs.cd[0,1],outwcs.wcs.cd[1,1]) * 180./np.pi
    else:
        outwcs.pscale = sqrt(outwcs.wcs.cd[0,0]**2 + outwcs.wcs.cd[1,0]**2)*3600.
        outwcs.orientat = arctan2(outwcs.wcs.cd[0,1],outwcs.wcs.cd[1,1]) * 180./np.pi
    
    tanpix = outwcs.wcs.s2p(fra_dec, 1)['pixcrd']
    
    outwcs.naxis1 = int(np.ceil(tanpix[:,0].max() - tanpix[:,0].min()))
    outwcs.naxis2 = int(np.ceil(tanpix[:,1].max() - tanpix[:,1].min()))
    crpix = np.array([outwcs.naxis1/2., outwcs.naxis2/2.])
    outwcs.wcs.crpix = crpix
    
    tanpix = outwcs.wcs.s2p(fra_dec, 1)['pixcrd']
    newcrpix = np.array([crpix[0]+np.ceil(tanpix[:,0].min()), crpix[1]+np.ceil(tanpix[:,1].min())])
    newcrval = outwcs.wcs.p2s([newcrpix], 1)['world'][0]
    outwcs.wcs.crval = newcrval
    return outwcs

def  undistortWCS(wcsobj):
    assert isinstance(wcsobj, pywcs.WCS)
    import coeff_converter
    
    cx, cy = coeff_converter.sip2idc(wcsobj)
    crpix1 = wcsobj.wcs.crpix[0]
    crpix2 = wcsobj.wcs.crpix[1]
    xy = np.array([(crpix1,crpix2),(crpix1+1.,crpix2),(crpix1,crpix2+1.)],dtype=np.double)
    offsets = np.array([wcsobj.ltv1, wcsobj.ltv2])
    px = xy + offsets
    #order = wcsobj.sip.a_order
    pscale = wcsobj.idcscale
    #pixref = np.array([wcsobj.sip.SIPREF1, wcsobj.sip.SIPREF2])
    
    tan_pix = apply_idc(px, cx, cy, wcsobj.wcs.crpix, pscale, order=1)
    xc = tan_pix[:,0]
    yc = tan_pix[:,1]
    am = xc[1] - xc[0]
    bm = xc[2] - xc[0]
    cm = yc[1] - yc[0]
    dm = yc[2] - yc[0]
    cd_mat = np.array([[am,bm],[cm,dm]],dtype=np.double)

    # Check the determinant for singularity
    _det = (am * dm) - (bm * cm)
    if ( _det == 0.0):
        print 'Singular matrix in updateWCS, aborting ...'
        return
    #lin_wcsobj = wcsutil.HSTWCS(instrument=wcsobj.instrument)
    lin_wcsobj = pywcs.WCS()    #instrument=wcsobj.instrument)
    cd_inv = np.linalg.inv(cd_mat)
    lin_wcsobj.wcs.cd = np.dot(wcsobj.wcs.cd, cd_inv)
    
    lin_wcsobj.orientat = arctan2(lin_wcsobj.wcs.cd[0,1],lin_wcsobj.wcs.cd[1,1]) * 180./np.pi
    lin_wcsobj.pscale = sqrt(lin_wcsobj.wcs.cd[0,0]**2 + lin_wcsobj.wcs.cd[1,0]**2)*3600.
    
    return lin_wcsobj

def apply_idc(pixpos, cx, cy, pixref, pscale= None, order=None):
    #pixpos must be already corrected for ltv1/2
    if cx == None:
        return pixpos  

    if order is None:
        print 'Unknown order of distortion model \n'
        return pixpos
    if pscale is None:
        print 'Unknown model plate scale\n'
        return pixpos
        
    # Apply in the same way that 'drizzle' would...
    _cx = cx/pscale
    _cy = cy/ pscale
    _p = pixpos
    
    # Do NOT include any zero-point terms in CX,CY here
    # as they should not be scaled by plate-scale like rest
    # of coeffs...  This makes the computations consistent
    # with 'drizzle'.  WJH 17-Feb-2004
    _cx[0,0] = 0.
    _cy[0,0] = 0.
    
    dxy = _p - pixref
    # Apply coefficients from distortion model here...
    
    c = _p * 0. 
    for i in range(order+1):
        for j in range(i+1):
            c[:,0] = c[:,0] + _cx[i][j] * pow(dxy[:,0],j) * pow(dxy[:,1],(i-j))
            c[:,1] = c[:,1] + _cy[i][j] * pow(dxy[:,0],j) * pow(dxy[:,1],(i-j))

    return  c


