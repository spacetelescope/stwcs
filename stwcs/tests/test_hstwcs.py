from astropy.io import fits
from stwcs.wcsutil import hstwcs, HSTWCS


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
