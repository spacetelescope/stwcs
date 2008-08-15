import wcsutil
from  wcsutil import HSTWCS
from pytools import fileutil, parseinput
import pyfits

#to avoid relative imports
import mappings
from mappings import basic_wcs
import distortion
import pywcs

def restoreWCS(fnames):
    """
    Given a list of fits file names restore the original basic WCS kw
    and write out the files. This overwrites the original files.
    """
    files = parseinput.parseinput(fnames)[0]
    for f in files:
        isfits, ftype = fileutil.isFits(f)
        if not isfits or (isfits and ftype == 'waiver'):
            print "RestoreWCS works only on multiextension fits files."
            return
        else:
            fobj = pyfits.open(f, mode='update')
            hdr0 = fobj[0].header
            for ext in fobj:
                try:
                    extname = ext.header['EXTNAME'].lower()
                except KeyError:
                    continue
                if extname in ['sci', 'err', 'sdq']:
                    hdr = ext.header
                    owcs = HSTWCS(hdr0, hdr)
                    #Changed restore to update the attributes ???
                    # this may need to be changed
                    backup = owcs.restore()
                    if not backup:
                        print '\No archived keywords found.\n'
                        continue
                    else:
                        for kw in basic_wcs:
                            nkw = ('O'+kw)[:7]
                            if backup.has_key(kw):
                                hdr.update(kw, hdr[nkw])
            fobj.close()
            
                    
