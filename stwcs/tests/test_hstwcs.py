import os

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from stwcs.wcsutil import HSTWCS, hstwcs


def test_radesys():
    phdr = fits.PrimaryHDU().header
    phdr['refframe'] = 'icrs'
    assert hstwcs.determine_refframe(phdr) == 'ICRS'
    phdr['refframe'] = 'gsc1'
    assert hstwcs.determine_refframe(phdr) == 'FK5'
    phdr['refframe'] = 'other'
    assert hstwcs.determine_refframe(phdr) is None
    phdr['refframe'] = ' '
    assert hstwcs.determine_refframe(phdr) is None


def test_prihdu_with_extver_no_extname():
    hdulist = fits.HDUList([
        fits.PrimaryHDU(header=fits.Header([('extver', 7)])),
        fits.ImageHDU(header=fits.Header([('time', 6)]))
    ])
    extname = HSTWCS(hdulist).extname
    assert extname == ('PRIMARY', 7)
    assert hdulist[extname] is hdulist[0]


def test_world_to_pixel():
    """
    Test scalar inputs for world to pixel conversions.
    """
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    fitsfile = os.path.join(data_dir, 'simple.fits')
    with fits.open(fitsfile) as hdulist:
        wcs = HSTWCS(hdulist, ext=0)

    ra = 5.5
    dec = -72.0
    sc_in = SkyCoord(ra=ra, dec=dec, unit='deg')
    x1, y1 = wcs.world_to_pixel(sc_in)
    x2, y2 = wcs.all_world2pix(ra, dec, 0)
    x3, y3 = wcs.wcs_world2pix(ra, dec, 0)
    for val in (x1, y1, x2, y2, x3, y3):
        assert isinstance(val, np.ndarray)
        assert val.shape == ()
        assert np.ndim(val) == 0
