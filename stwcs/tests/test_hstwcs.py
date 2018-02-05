from astropy.io import fits
from ..wcsutil import hstwcs


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
