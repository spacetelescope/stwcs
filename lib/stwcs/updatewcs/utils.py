from __future__ import division # confidence high
import os

from stsci.tools import fileutil

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

def extract_rootname(kwvalue,suffix=""):
    """ Returns the rootname from a full reference filename

        If a non-valid value (any of ['','N/A','NONE','INDEF',None]) is input,
            simply return a string value of 'NONE'
            
        This function will also replace any 'suffix' specified with a blank.
    """
    # check to see whether a valid kwvalue has been provided as input
    if kwvalue.strip() in ['','N/A','NONE','INDEF',None]:
        return 'NONE' # no valid value, so return 'NONE'
    
    # for a valid kwvalue, parse out the rootname
    # strip off any environment variable from input filename, if any are given
    if '$' in kwvalue:
        fullval = kwvalue[kwvalue.find('$')+1:]
    else:
        fullval = kwvalue
    # Extract filename without path from kwvalue
    fname = os.path.basename(fullval).strip()
    
    # Now, rip out just the rootname from the full filename
    rootname = fileutil.buildNewRootname(fname)

    # Now, remove any known suffix from rootname
    rootname = rootname.replace(suffix,'')
    return rootname.strip()

def construct_distname(fobj,wcsobj):
    """ 
    This function constructs the value for the keyword 'DISTNAME'. 
    It relies on the reference files specified by the keywords 'IDCTAB',
    'NPOLFILE', and 'D2IMFILE'.  
    
    The final constructed value will be of the form:
        <idctab rootname>-<npolfile rootname>-<d2imfile rootname>
    and have a value of 'NONE' if no reference files are specified.
    """
    idcname = extract_rootname(fobj[0].header.get('IDCTAB', "NONE"),
                                suffix='_idc')
    if idcname is None and wcsobj.sip is not None:
        idcname = 'UNKNOWN'
        
    npolname = build_npolname(fobj)
    if npolname is None and wcsobj.cpdis1 is not None:
        npolname = 'UNKNOWN'
    
    d2imname = build_d2imname(fobj)
    if d2imname is None and wcsobj.det2im is not None:
        d2imname = 'UNKNOWN'
    
    sipname = build_sipname(fobj)

    distname = build_distname(sipname,npolname,d2imname)
    return {'DISTNAME':distname,'SIPNAME':sipname}

def build_distname(sipname,npolname,d2imname):
    """ 
    Core function to build DISTNAME keyword value without the HSTWCS input.
    """

    distname = sipname.strip()
    if npolname != 'NONE' or d2imname != 'NONE':
        if d2imname != 'NONE':
            distname+= '-'+npolname.strip() + '-'+d2imname.strip()
        else:
            distname+='-'+npolname.strip()

    return distname

def build_default_wcsname(idctab):

    idcname = extract_rootname(idctab,suffix='_idc')
    wcsname = 'IDC_' + idcname
    return wcsname

def build_sipname(fobj):
    try:
        sipname = fobj[0].header["SIPNAME"]
    except KeyError:
        try:
            idcname = extract_rootname(fobj[0].header["IDCTAB"],suffix='_idc')
            sipname = fobj[0].header['rootname']+'_'+ idcname
        except KeyError:
            if 'A_ORDER' in fobj[1].header or 'B_ORDER' in fobj[1].header:
                sipname = 'UNKNOWN'
            else:
                sipname = 'NOMODEL'
    return sipname

def build_npolname(fobj):
    try:
        npolfile = extract_rootname(fobj[0].header["NPOLFILE"],suffix='_npl')
    except KeyError:
        if fileutil.countExtn(fobj, 'WCSDVARR'):
            npolfile = 'UNKNOWN'
        else:
            npolfile = 'NOMODEL'
    return npolfile

def build_d2imname(fobj):
    try:
        d2imfile = extract_rootname(fobj[0].header["D2IMFILE"],suffix='_d2i')
    except KeyError:
        if fileutil.countExtn(fobj, 'D2IMARR'):
            d2imfile = 'UNKNOWN'
        else:
            d2imfile = 'NOMODEL'
    return d2imfile