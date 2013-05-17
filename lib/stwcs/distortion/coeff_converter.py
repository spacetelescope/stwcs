from __future__ import division # confidence high

import numpy as np
from astropy.io import fits
from astropy import wcs

def sip2idc(hwcs):
    """
    Converts SIP style coefficients to IDCTAB coefficients.

    :Parameters:
    `wcs`: fits.Header or pywcs.WCS object
    """
    if isinstance(hwcs,fits.Header):
        ocx10 = hwcs.get('OCX10', None)
        ocx11 = hwcs.get('OCX11', None)
        ocy10 = hwcs.get('OCY10', None)
        ocy11 = hwcs.get('OCY11', None)
        order = hwcs.get('A_ORDER', None)
        sipa, sipb = _read_sip_kw(hwcs)
        if None in [ocx10, ocx11, ocy10, ocy11, sipa, sipb]:
            print 'Cannot convert SIP to IDC coefficients.\n'
            return None, None
    elif isinstance(hwcs,wcs.WCS):
        try:
            ocx10 = hwcs.ocx10
            ocx11 = hwcs.ocx11
            ocy10 = hwcs.ocy10
            ocy11 = hwcs.ocy11
        except AttributeError:
            print 'First order IDCTAB coefficients are not available.\n'
            print 'Cannot convert SIP to IDC coefficients.\n'
            return None, None
        try:
            sipa = hwcs.sip.a
            sipb = hwcs.sip.b
        except AttributeError:
            print 'SIP coefficients are not available.'
            print 'Cannot convert SIP to IDC coefficients.\n'
            return None, None
        try:
            order = hwcs.sip.a_order
        except AttributeError:
            print 'SIP model order unknown, exiting ...\n'
            return None, None

    else:
        print 'Input to sip2idc must be a FITS header or a wcsutil.HSTWCS object\n'
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

def _read_sip_kw(header):
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


