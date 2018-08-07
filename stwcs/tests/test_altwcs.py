import shutil
import os
from astropy.io import fits as pyfits
from ..wcsutil import altwcs
from .. import updatewcs
from .. import wcsutil
import numpy as np
from numpy.testing import utils
import pytest

from . import data
data_path = os.path.split(os.path.abspath(data.__file__))[0]

os.environ['ASTROMETRY_STEP_CONTROL'] = 'Off'

def get_filepath(filename, directory=data_path):
    return os.path.join(directory, filename)


def compare_wcs(w1, w2, exclude_keywords=None):
    """
    Compare two WCSs.

    Parameters
    ----------
    w1, w2 : `astropy.wcs.WCS` objects
    exclude_keywords : list
        List of keywords to excude from comparison.
    """
    exclude_ctype = False
    keywords = ['crval', 'crpix', 'cd']
    if exclude_keywords is not None:
        exclude_keywords = [kw.lower() for kw in exclude_keywords]
        if 'ctype' in exclude_keywords:
            exclude_ctype = True
            exclude_keywords.remove('ctype')
        for kw in exclude_keywords:
            keywords.remove(kw)
    for kw in keywords:
        kw1 = getattr(w1.wcs, kw)
        kw2 = getattr(w2.wcs, kw)
        utils.assert_allclose(kw1, kw2, 1e-10)
    #utils.assert_allclose(w1.wcs.crpix, w2.wcs.crpix, 1e-10)
    #utils.assert_allclose(w1.wcs.cd, w2.wcs.cd, 1e-10)
    if not exclude_ctype:
        utils.assert_array_equal(np.array(w1.wcs.ctype), np.array(w2.wcs.ctype))

class TestAltWCS(object):

    def setup_class(self):
        acs_orig_file = get_filepath('j94f05bgq_flt.fits')
        simple_orig_file = get_filepath('simple.fits')
        current_dir = os.path.abspath(os.path.curdir)
        simple_file = get_filepath('simple.fits', current_dir)
        acs_file = get_filepath('j94f05bgq_flt.fits', current_dir)

        try:
            os.remove(acs_file)
            os.remove(simple_file)
        except OSError:
            pass
        idctab = get_filepath('postsm4_idc.fits')
        npol_file = get_filepath('qbu16424j_npl.fits')
        d2imfile = get_filepath('new_wfc_d2i.fits ')

        shutil.copyfile(acs_orig_file, acs_file)
        shutil.copyfile(simple_orig_file, simple_file)
        pyfits.setval(acs_file, ext=0, keyword="IDCTAB", value=idctab)
        pyfits.setval(acs_file, ext=0, keyword="NPOLFILE", value=npol_file)
        pyfits.setval(acs_file, ext=0, keyword="D2IMFILE", value=d2imfile)

        updatewcs.updatewcs(acs_file)
        self.acs_file = acs_file
        self.simplefits = simple_file
        self.ww = wcsutil.HSTWCS(self.acs_file, ext=1)

    def test_archive(self):
        altwcs.archiveWCS(self.acs_file, ext=1, wcskey='Z', wcsname='ZTEST', reusekey=False)
        w1 = wcsutil.HSTWCS(self.acs_file, ext=1)
        w1z = wcsutil.HSTWCS(self.acs_file, ext=1, wcskey='Z')
        compare_wcs(w1, w1z)

    def test_archive_clobber(self):
        altwcs.archiveWCS(self.acs_file, ext=1, wcskey='Z', wcsname='ZTEST', reusekey=True)
        w1 = wcsutil.HSTWCS(self.acs_file, ext=1)
        w1z = wcsutil.HSTWCS(self.acs_file, ext=1, wcskey='Z')
        compare_wcs(w1, w1z)

    def test_restore_wcs(self):
        # test restore on a file
        altwcs.restoreWCS(self.acs_file, ext=1, wcskey='O')
        w1o = wcsutil.HSTWCS(self.acs_file, ext=1, wcskey='O')
        w1 = wcsutil.HSTWCS(self.acs_file, ext=1)
        compare_wcs(w1, w1o, exclude_keywords=['ctype'])

    def test_restore_wcs_mem(self):
        # test restore on an HDUList object
        altwcs.archiveWCS(self.acs_file, ext=[('SCI', 1), ('SCI', 2)], wcskey='T')
        pyfits.setval(self.acs_file, ext=('SCI', 1), keyword='CRVAL1', value=1)
        pyfits.setval(self.acs_file, ext=('SCI', 2), keyword='CRVAL1', value=1)
        f = pyfits.open(self.acs_file, mode='update')
        altwcs.restoreWCS(f, ext=1, wcskey='T')
        f.close()
        w1o = wcsutil.HSTWCS(self.acs_file, ext=1, wcskey='T')
        w1 = wcsutil.HSTWCS(self.acs_file, ext=1)
        compare_wcs(w1, w1o)

    def test_restore_simple(self):
        # test restore on simple fits format
        altwcs.archiveWCS(self.simplefits, ext=0, wcskey='R')
        pyfits.setval(self.simplefits, ext=0, keyword='CRVAL1R', value=1)
        altwcs.restoreWCS(self.simplefits, ext=0, wcskey='R')
        wo = wcsutil.HSTWCS(self.simplefits, ext=0, wcskey='R')
        ws = wcsutil.HSTWCS(self.simplefits, ext=0)
        compare_wcs(ws, wo)

    def test_restore_wcs_from_to(self):
        # test restore from ... to ...
        pyfits.setval(self.acs_file, ext=('SCI', 1), keyword='CRVAL1', value=1)
        pyfits.setval(self.acs_file, ext=('SCI', 2), keyword='CRVAL1', value=1)
        f = pyfits.open(self.acs_file, mode='update')
        altwcs.restore_from_to(f, fromext='SCI', toext=['SCI', 'ERR', 'DQ'],
                               wcskey='T')
        f.close()
        w1o = wcsutil.HSTWCS(self.acs_file, ext=('SCI', 1), wcskey='T')
        w1 = wcsutil.HSTWCS(self.acs_file, ext=('SCI', 1))
        compare_wcs(w1, w1o)
        w2 = wcsutil.HSTWCS(self.acs_file, ext=('ERR', 1))
        compare_wcs(w2, w1o, exclude_keywords=['ctype'])
        w3 = wcsutil.HSTWCS(self.acs_file, ext=('DQ', 1))
        compare_wcs(w3, w1o, exclude_keywords=['ctype'])
        w4o = wcsutil.HSTWCS(self.acs_file, ext=4, wcskey='T')
        w4 = wcsutil.HSTWCS(self.acs_file, ext=('SCI', 2))
        compare_wcs(w4, w4o)
        w5 = wcsutil.HSTWCS(self.acs_file, ext=('ERR', 2))
        compare_wcs(w5, w4o, exclude_keywords=['ctype'])
        w6 = wcsutil.HSTWCS(self.acs_file, ext=('DQ', 2))
        compare_wcs(w3, w1o, exclude_keywords=['ctype'])

    def test_delete_wcs(self):
        altwcs.deleteWCS(self.acs_file, ext=1, wcskey='Z')
        with pytest.raises(KeyError):
            wcsutil.HSTWCS(self.acs_file, ext=1, wcskey='Z')

    def test_pars_file_mode1(self):
        assert(not altwcs._parpasscheck(self.acs_file, ext=1, wcskey='Z'))

    def test_pars_file_mode2(self):
        f = pyfits.open(self.acs_file)
        assert(not altwcs._parpasscheck(f, ext=1, wcskey='Z'))
        f.close()

    def test_pars_ext(self):
        f = pyfits.open(self.acs_file, mode='update')
        assert(altwcs._parpasscheck(f, ext=1, wcskey='Z'))
        assert(altwcs._parpasscheck(f, ext=[('sci', 1), ('sci', 2)], wcskey='Z'))
        assert(altwcs._parpasscheck(f, ext=('sci', 1), wcskey='Z'))
        f.close()

    def test_pars_wcskey_not1char(self):
        f = pyfits.open(self.acs_file, mode='update')
        assert(not altwcs._parpasscheck(f, ext=1, wcskey='ZZ'))
        f.close()

    def test_pars_wcskey(self):
        f = pyfits.open(self.acs_file, mode='update')
        assert(altwcs._parpasscheck(f, ext=1, wcskey=' '))
        #assert(not altwcs._parpasscheck(f, ext=1, wcskey=' ', reusekey=False))
        #assert(altwcs._parpasscheck(f, ext=1, wcskey='O'))
        #assert(not altwcs._parpasscheck(f, ext=1, wcskey='O', reusekey=False))
        f.close()
