from __future__ import division # confidence high
import os

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

def extract_rootname(kwvalue):
    """ Returns the rootname from a full reference filename

        If a non-valid value (any of ['','N/A','NONE','INDEF',None]) is input,
            simply return a string value of 'NONE'
            
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
    if '_' in fname:
        rootname = fname[:fname.rfind('_')]
    elif '.fit' in fname: # on some systems, only .fit is used instead of .fits
        rootname = fname[:fname.rfind('.fit')]
    elif '.' in fname: # account for non-standard file extensions
        rootname = fname[:fname.rfind('.')]
    else:
        rootname = fname
        
    return rootname

def construct_distname(fobj,wcsobj):
    """ 
    This function constructs the value for the keyword 'DISTNAME'. 
    It relies on the reference files specified by the keywords 'IDCTAB',
    'NPOLFILE', and 'D2IMFILE'.  
    
    The final constructed value will be of the form:
        <idctab rootname>-<npolfile rootname>-<d2imfile rootname>
    and have a value of 'NONE' if no reference files are specified.
    """
    idcname = extract_rootname(fobj[0].header.get('IDCTAB', "NONE"))
    if idcname is None and wcsobj.sip is not None:
        idcname = 'UNKNOWN'
        
    npolname = extract_rootname(fobj[0].header.get('NPOLFILE', "NONE"))
    if npolname is None and wcsobj.cpdis1 is not None:
        npolname = 'UNKNOWN'
    
    d2imname = extract_rootname(fobj[0].header.get('D2IMFILE', "NONE"))
    if d2imname is None and wcsobj.det2im is not None:
        d2imname = 'UNKNOWN'
    
    sipname = '%s_%s'%(fobj[0].header.get('rootname',""),idcname)
    distname = sipname.strip()
    if npolname != 'NONE' or d2imname != 'NONE':
        if d2imname != 'NONE':
            distname+= '-'+npolname.strip() + '-'+d2imname.strip()
        else:
            distname+='-'+npolname.strip()

    return {'DISTNAME':distname,'SIPNAME':sipname}