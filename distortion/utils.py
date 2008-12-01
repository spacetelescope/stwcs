import numpy as np
import pywcs
import pyfits
from hstwcs import wcsutil
from numpy import sqrt, arctan2

def output_wcs(list_of_wcsobj, ref_wcs=None, outwcs=None):
    fra_dec = np.vstack([w.footprint for w in list_of_wcsobj])
    """
    ra_min = np.array(fra_dec[:,0]).min()
    dec_min = np.array(fra_dec[:,1]).min()
    ra_max = np.array(fra_dec[:,0]).max()
    dec_max = np.array(fra_dec[:,1]).max()
                   
    output_footprint=np.zeros(shape=(4,2),dtype=np.float64)
    output_footprint[0,0]=ra_min
    output_footprint[0,1]=dec_min
    output_footprint[1,0]=ra_min
    output_footprint[1,1]=dec_max
    output_footprint[2,0]=ra_max
    output_footprint[2,1]=dec_max
    output_footprint[3,0]=ra_max
    output_footprint[3,1]=dec_min
    """
    if outwcs is None:
        if ref_wcs == None:
            ref_wcs = list_of_wcsobj[0]
            
        outwcs = undistortWCS(ref_wcs)
        outwcs.wcs.crpix = ref_wcs.wcs.crpix
        outwcs.wcs.crval = ref_wcs.wcs.crval
        outwcs.pscale = sqrt(outwcs.wcs.cd[0,0]**2 + outwcs.wcs.cd[1,0]**2)*3600.
        outwcs.orientat = arctan2(outwcs.wcs.cd[0,1],outwcs.wcs.cd[1,1]) * 180./np.pi
    else:
        outwcs.pscale = sqrt(outwcs.wcs.cd[0,0]**2 + outwcs.wcs.cd[1,0]**2)*3600.
        outwcs.orientat = arctan2(outwcs.wcs.cd[0,1],outwcs.wcs.cd[1,1]) * 180./np.pi
    
    #out_px = outwcs.wcs.s2p_fits(output_footprint)['pixcrd']
    out_px = outwcs.wcs.s2p_fits(fra_dec)['pixcrd']
    outwcs.naxis1 = int(np.ceil(out_px[:,0].max()) - np.floor(out_px[:,0].min()))
    outwcs.naxis2 = int(np.ceil(out_px[:,1].max()) - np.floor(out_px[:,1].min()))
    outwcs.recenter()
    
    return outwcs

def  undistortWCS(wcsobj):
    assert isinstance(wcsobj, pywcs.WCS)
    #if wcsobj.idcmodel == None:
    #    return

    # get the sip coefficients and the first order IDC coeffs
    # reconstruct the idc model
    # apply the idc model
    
    from hstwcs.distortion import coeff_converter
    
    cx, cy = coeff_converter.sip2idc(wcsobj)
    crpix1 = wcsobj.wcs.crpix[0]
    crpix2 = wcsobj.wcs.crpix[1]
    xy = np.array([(crpix1,crpix2),(crpix1+1.,crpix2),(crpix1,crpix2+1.)],dtype=np.double)
    offsets = np.array([wcsobj.ltv1, wcsobj.ltv2])
    px = xy + offsets
    order = wcsobj.sip.a_order
    pscale = wcsobj.idcscale
    #pixref = np.array([wcsobj.sip.SIPREF1, wcsobj.sip.SIPREF2])
    
    tan_pix = apply_idc(px, cx, cy, wcsobj.wcs.crpix, pscale, order)
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
    print 'inv_cd', cd_inv
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

def sip2idc(wcs):
    """
    Converts SIP style coefficients to IDCTAB coefficients.
    
    :Parameters:
    `wcs`: pyfits.Header or pywcs.WCS object
    """
    if isinstance(wcs,pyfits.Header):
        ocx10 = wcs.get('OCX10', None)
        ocx11 = wcs.get('OCX11', None)
        ocy10 = wcs.get('OCY10', None)
        ocy11 = wcs.get('OCY11', None)
        order = hdr.get('A_ORDER', None)
        sipa, sipb = read_sip_kw(header)
        if sipa == None or sipb == None:
            print 'SIP coefficients are not available.\n'
            print 'Cannot convert SIP to IDC coefficients.\n'
            return
    elif isinstance(wcs,pywcs.WCS):
        try:
            ocx10 = wcs.ocx10
            ocx11 = wcs.ocx11
            ocy10 = wcs.ocy10
            ocy11 = wcs.ocy11
        except AttributeError:
            print 'First order IDCTAB coefficients are not available.\n'
            print 'Cannot convert SIP to IDC coefficients.\n'
            return
        try:
            sipa = wcs.sip.a
            sipb = wcs.sip.b
        except AttributeError:
            print 'SIP coefficients are not available.\n'
            print 'Cannot convert SIP to IDC coefficients.\n'
            return
    else:
        print 'Input to sip2idc must be a PyFITS header or a wcsutil.HSTWCS object\n'
        return
    
    try:
        order = wcs.sip.a_order
    except AttributeError:
        print 'SIP model order unknown, exiting ...\n'
        return

    if None in [ocx10, ocx11, ocy10, ocy11]:
        print 'First order IDC coefficients not found, exiting ...\n'
        return
    idc_coeff = np.array([[ocx11, ocx10], [ocy11, ocy10]])
    cx = np.zeros((order+1,order+1), dtype=np.double)
    cy = np.zeros((order+1,order+1), dtype=np.double)
    for n in range(order+1):
        for m in range(order+1):
            if n >= m and n>=2:
                sipval = np.array([[sipa[m,n-m]],[sipb[m,n-m]]])
                idcval = np.dot(idc_coeff, sipval)
                cx[n,m] = idcval[0]
                cy[n,m] = idcval[1]
    
    cx[1,0] = ocx10
    cx[1,1] = ocx11
    cy[1,0] = ocy10
    cy[1,1] = ocy11
                
    return cx, cy

def read_sip_kw(header):
    """
    Reads SIP header keywords and returns an array of coefficients.

    If no SIP header keywords are found, None is returned.
    """
    if header.has_key("A_ORDER"):
        if not header.has_key("B_ORDER"):
            raise ValueError(
                "A_ORDER provided without corresponding B_ORDER "
                "keyword for SIP distortion")

        m = int(header["A_ORDER"])
        a = np.zeros((m+1, m+1), np.double)
        for i in range(m+1):
            for j in range(m-i+1):
                a[i, j] = header.get("A_%d_%d" % (i, j), 0.0)

        m = int(header["B_ORDER"])
        b = np.zeros((m+1, m+1), np.double)
        for i in range(m+1):
            for j in range(m-i+1):
                b[i, j] = header.get("B_%d_%d" % (i, j), 0.0)
    elif header.has_key("B_ORDER"):
        raise ValueError(
            "B_ORDER provided without corresponding A_ORDER "
            "keyword for SIP distortion")
    else:
        a = None
        b = None
            
    return a , b

