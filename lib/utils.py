from __future__ import division # confidence high

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

def getBinning(fobj, extver=1):
    # Return the binning factor
    binned = 1
    if fobj[0].header['INSTRUME'] == 'WFPC2':
        mode = fobj[0].header.get('MODE', "")
        if mode == 'AREA': binned = 2
    else:
        binned = fobj['SCI', extver].header.get('BINAXIS',1)
    return binned

