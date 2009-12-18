from pytools import parseinput, fileutil
import pyfits
from wcsutil.mappings import basic_wcs

def restoreWCS(fnames):
    """
    Given a list of fits file names restore the original basic WCS kw
    and write out the files. This overwrites the original files.
    
    Affected keywords:
    'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CRVAL1','CRVAL2','CTYPE1', 'CTYPE2',
    'CRPIX1', 'CRPIX2', 'CTYPE1', 'CTYPE2', 'ORIENTAT'
    """
    files = parseinput.parseinput(fnames)[0]
    for f in files:
        isfits, ftype = fileutil.isFits(f)
        if not isfits or (isfits and ftype == 'waiver'):
            print "RestoreWCS works only with true fits files."
            return
        else:
            fobj = pyfits.open(f, mode='update')
            for ext in fobj:
                try:
                    extname = ext.header['EXTNAME'].lower()
                except KeyError:
                    continue
                if extname in ['sci', 'err', 'sdq']:
                    hdr = ext.header
                    backup = get_archive(hdr)
                    if not backup:
                        #print 'No archived keywords found.\n'
                        continue
                    else:
                        for kw in basic_wcs:
                            nkw = ('O'+kw)[:7]
                            if backup.has_key(kw):
                                hdr.update(kw, hdr[nkw])
                    tddalpha = hdr.get('TDDALPHA', None)
                    tddbeta = hdr.get('TDDBETA', None)
                    if tddalpha or tddbeta:
                        hdr.update('TDDALPHA', 0.0)
                        hdr.update('TDDBETA', 0.0)
            fobj.close()

def get_archive(header):
    """
    Returns a dictionary with the archived basic WCS keywords.
    """    

    backup = {}
    for k in basic_wcs:
        try:
            nkw = ('O'+k)[:7]
            backup[k] = header[nkw]
        except KeyError:
            pass
    return backup
    
def write_archive(header):
    """
    Archives original WCS kw before recalculating them.
    """
    backup_kw = get_archive(header)
    if backup_kw != {}:
        #print "Archive already exists\n"
        return
    else:
        cmt = 'archived value'
        for kw in basic_wcs:
            nkw = 'O'+kw
            try:
                header.update(nkw[:8], header[kw], comment=cmt)
            except KeyError:  #ORIENTAT is not always present in headers
                pass
        header.update('WCSCDATE', fileutil.getLTime(), comment='Local time the WCS kw were archived')
    

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
