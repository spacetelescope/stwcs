from astropy.io import fits
from stwcs.wcsutil import hstwcs, HSTWCS
import numpy as np


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


def test_pc2cd():
    e = np.exp(1)
    hdulist = fits.HDUList([
        fits.PrimaryHDU(header=fits.Header([('extver', 1)])),
        fits.ImageHDU(header=fits.Header([
            ('pc1_1', 0.1),
            ('pc1_2', 0.08),
            ('pc2_1', -0.05),
            ('pc2_2', 0.7),
            ('cdelt1', e),
            ('cdelt2', np.pi),
            ('crpix1', 512.0),
            ('crpix2', 512.0),
            ('crval1', 150.0),
            ('crval2', 2.0),
            ('ctype1', 'RA---TAN'),
            ('ctype2', 'DEC--TAN'),
        ]))
    ])
    wcs = HSTWCS(hdulist, ext=1)
    assert wcs.wcs.has_cd()
    cd = np.array([[0.1 * e, 0.08 * e], [-0.05 * np.pi, 0.7 * np.pi]])
    assert np.allclose(wcs.wcs.cd, cd)
