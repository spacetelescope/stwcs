#!/usr/bin/env python

from pydrizzle import process_input
import pywcs, pyfits, numpy as np
from stwcs.distortion import utils
from pytools import fileutil as fu
from stwcs import wcsutil


def align_wcs(input, shiftfile=None, writeasn=False, asnname=None):
    """
    Purpose
    =======
    Correct the WCS of a fits file, so that multidrizzle aligns the images.
    To view the format of the shift file:
    >>>from pytools.asnutil import ShiftFile  
    >>>print ShiftFile.__doc__

    Example
    =======
    >>>import alignwcs
    >>>alignwcs.align_wcs('*flt.fits', shiftfile='shifts.txt')
    It works also on the command line:
    %./alignwcs.py *flt.fits shifts.txt
    
    Dependencies 
    ============
    `pytools`
    `pyfits`
    `pywcs`
    `numpy`
    `pydrizzle`

    :Parameters:
    `input`: a python list of file names or a string (wild card characters allowed)
             input files may be in fits, geis or waiver fits format
    `shiftfile`: a shift file, as the output of tweakshifts
    `writeasn`: boolean
                defines whether to write out the asn table updated with shift information
    `asnname`: string
                name for the output asn table
    """
    if shiftfile == None:
        print 'A shift file is required but was not provided\n'
        return
    
    asndict,b,c = process_input.process_input(input, shiftfile=shiftfile, updatewcs=False)
    if writeasn:
        if asnname == None:
            asnname = c.split('.fits')[0] + '_shifts' + '.fits'
        a.write(asnname)
    outwcs = get_outwcs(asndict)
    apply_shifts(asndict, outwcs)

def apply_shifts(asndict, outwcs):
    for mem in asndict['members']:
        filename = fu.buildRootname(mem)
        xsh = asndict['members'][mem]['xshift']
        ysh = asndict['members'][mem]['yshift']
        shifts = np.array([[xsh, ysh]])
        f = pyfits.open(filename)
        for extn in f:
            if extn.header.has_key('extname') and extn.header['extname'].lower() == 'sci':
                extver = extn.header['extver']
                owcs = pywcs.WCS(extn.header, f)
                crvals = np.array([owcs.wcs.crval])
                px = outwcs.wcs.s2p(crvals, 1)['pixcrd'] + shifts
                ncrvals = outwcs.all_pix2sky(px, 1)
                pyfits.setval(filename, 'CRVAL1', value=ncrvals[0,0], ext=extver)
                pyfits.setval(filename, 'CRVAL2', value=ncrvals[0,1], ext=extver)
                print 'Updated %s with shifts ' % filename, shifts
        f.close()


def get_outwcs(asndict):
    wcslist = []
    for mem in asndict['members']:
        filename = fu.buildRootname(mem)
        f = pyfits.open(filename)
        hdr0 = f[0].header
        for extn in f:
            if extn.header.has_key('extname') and extn.header['extname'].lower() == 'sci':
                owcs = wcsutil.HSTWCS(hdr0, extn.header,f)
                wcslist.append(owcs)
        f.close()
    outwcs = utils.output_wcs(wcslist)
    return outwcs

if __name__ == '__main__':
    import sys
    args = sys.argv[1:]
    input = args[:-1]
    shifts = args[-1]
    
    align_wcs(input, shiftfile=shifts)
