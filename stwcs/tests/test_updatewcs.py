import shutil
import os

from astropy import wcs
from astropy.io import fits
from .. import updatewcs
from ..updatewcs import apply_corrections
from ..distortion import utils as dutils
from ..wcsutil import HSTWCS
import numpy as np
from numpy.testing import utils
import pytest


from . import data
data_path = os.path.split(os.path.abspath(data.__file__))[0]

os.environ['ASTROMETRY_STEP_CONTROL'] = 'Off'

def get_filepath(filename, directory=data_path):
    return os.path.join(directory, filename)


class TestStwcs(object):

    def setup_class(self):
        acs_orig_file = get_filepath('j94f05bgq_flt.fits')
        current_dir = os.path.abspath(os.path.curdir)
        self.ref_file = get_filepath('j94f05bgq_flt_r.fits', current_dir)
        self.acs_file = get_filepath('j94f05bgq_flt.fits', current_dir)

        try:
            os.remove(self.acs_file)
            os.remove(self.ref_file)
        except OSError:
            pass
        idctab = get_filepath('postsm4_idc.fits')
        npol_file = get_filepath('qbu16424j_npl.fits')
        d2imfile = get_filepath('new_wfc_d2i.fits ')

        shutil.copyfile(acs_orig_file, self.acs_file)

        fits.setval(self.acs_file, ext=0, keyword="IDCTAB", value=idctab)
        fits.setval(self.acs_file, ext=0, keyword="NPOLFILE", value=npol_file)
        fits.setval(self.acs_file, ext=0, keyword="D2IMFILE", value=d2imfile)
        updatewcs.updatewcs(self.acs_file)
        #self.ref_file = ref_file
        shutil.copyfile(self.acs_file, self.ref_file)

        self.w1 = HSTWCS(self.acs_file, ext=1)
        self.w4 = HSTWCS(self.acs_file, ext=4)
        self.w = wcs.WCS()

    def test_default(self):
        crval = np.array([0., 0.])
        crpix = np.array([0., 0.])
        cdelt = np.array([1., 1.])
        pc = np.array([[1., 0], [0., 1.]])
        ctype = np.array(['', ''])
        utils.assert_almost_equal(self.w.wcs.crval, crval)
        utils.assert_almost_equal(self.w.wcs.crpix, crpix)
        utils.assert_almost_equal(self.w.wcs.cdelt, cdelt)
        utils.assert_almost_equal(self.w.wcs.pc, pc)
        assert((self.w.wcs.ctype == np.array(['', ''])).all())

    def test_simple_sci1(self):
        """
        A simple sanity check that CRPIX corresponds to CRVAL within wcs
        """
        px1 = np.array([self.w1.wcs.crpix])
        rd1 = np.array([self.w1.wcs.crval])
        assert(((self.w1.all_pix2world(px1, 1) - rd1) < 5E-7).all())

    def test_simple_sci2(self):
        """
        A simple sanity check that CRPIX corresponds to CRVAL within wcs
        """
        px4 = np.array([self.w4.wcs.crpix])
        rd4 = np.array([self.w4.wcs.crval])
        assert(((self.w4.all_pix2world(px4, 1) - rd4) < 2E-6).all())

    def test_pipeline_sci1(self):
        """
        Internal consistency check of the wcs pipeline
        """
        px = np.array([[100, 125]])
        sky1 = self.w1.all_pix2world(px, 1)
        dpx1 = self.w1.det2im(px, 1)
        #fpx1 = dpx1 + (self.w1.sip_pix2foc(dpx1,1)-dpx1) + (self.w1.p4_pix2foc(dpx1,1)-dpx1)
        fpx1 = dpx1 + (self.w1.sip_pix2foc(dpx1, 1)-dpx1+self.w1.wcs.crpix) + \
            (self.w1.p4_pix2foc(dpx1, 1)-dpx1)
        pipelinepx1 = self.w1.wcs_pix2world(fpx1, 1)
        utils.assert_almost_equal(pipelinepx1, sky1)

    def test_pipeline_sci2(self):
        """
        Internal consistency check of the wcs pipeline
        """
        px = np.array([[100, 125]])
        sky4 = self.w4.all_pix2world(px, 1)
        dpx4 = self.w4.det2im(px, 1)
        fpx4 = dpx4 + (self.w4.sip_pix2foc(dpx4, 1)-dpx4 + self.w4.wcs.crpix) + \
            (self.w4.p4_pix2foc(dpx4, 1)-dpx4)
        pipelinepx4 = self.w4.wcs_pix2world(fpx4, 1)
        utils.assert_almost_equal(pipelinepx4, sky4)

    def test_outwcs(self):
        """
        Test the WCS of the output image
        """
        outwcs = dutils.output_wcs([self.w1, self.w4])

        #print('outwcs.wcs.crval = {0}'.format(outwcs.wcs.crval))
        utils.assert_allclose(
            outwcs.wcs.crval, np.array([5.65109952, -72.0674181]), atol=1e-7)

        utils.assert_almost_equal(outwcs.wcs.crpix, np.array([2107.0, 2118.5]))
        utils.assert_almost_equal(
            outwcs.wcs.cd,
            np.array([[1.2787045268089949e-05, 5.4215042082174661e-06],
                      [5.4215042082174661e-06, -1.2787045268089949e-05]]))
        assert(outwcs._naxis1 == 4214)
        assert(outwcs._naxis2 == 4237)

    def test_repeate(self):
        # make sure repeated runs of updatewcs do not change the WCS.
        updatewcs.updatewcs(self.acs_file)
        w1 = HSTWCS(self.acs_file, ext=('SCI', 1))
        w4 = HSTWCS(self.acs_file, ext=('SCI', 2))
        w1r = HSTWCS(self.ref_file, ext=('SCI', 1))
        w4r = HSTWCS(self.ref_file, ext=('SCI', 2))
        utils.assert_almost_equal(w1.wcs.crval, w1r.wcs.crval)
        utils.assert_almost_equal(w1.wcs.crpix, w1r.wcs.crpix)
        utils.assert_almost_equal(w1.wcs.cdelt, w1r.wcs.cdelt)
        utils.assert_almost_equal(w1.wcs.cd, w1r.wcs.cd)
        assert((np.array(w1.wcs.ctype) == np.array(w1r.wcs.ctype)).all())
        utils.assert_almost_equal(w1.sip.a, w1r.sip.a)
        utils.assert_almost_equal(w1.sip.b, w1r.sip.b)
        utils.assert_almost_equal(w4.wcs.crval, w4r.wcs.crval)
        utils.assert_almost_equal(w4.wcs.crpix, w4r.wcs.crpix)
        utils.assert_almost_equal(w4.wcs.cdelt, w4r.wcs.cdelt)
        utils.assert_almost_equal(w4.wcs.cd, w4r.wcs.cd)
        assert((np.array(self.w4.wcs.ctype) == np.array(w4r.wcs.ctype)).all())
        utils.assert_almost_equal(w4.sip.a, w4r.sip.a)
        utils.assert_almost_equal(w4.sip.b, w4r.sip.b)


def test_remove_npol_distortion():
    acs_orig_file = get_filepath('j94f05bgq_flt.fits')
    current_dir = os.path.abspath(os.path.curdir)
    acs_file = get_filepath('j94f05bgq_flt.fits', current_dir)
    idctab = get_filepath('postsm4_idc.fits')
    npol_file = get_filepath('qbu16424j_npl.fits')
    d2imfile = get_filepath('new_wfc_d2i.fits ')

    try:
        os.remove(acs_file)
    except OSError:
        pass

    shutil.copyfile(acs_orig_file, acs_file)
    fits.setval(acs_file, ext=0, keyword="IDCTAB", value=idctab)
    fits.setval(acs_file, ext=0, keyword="NPOLFILE", value=npol_file)
    fits.setval(acs_file, ext=0, keyword="D2IMFILE", value=d2imfile)

    updatewcs.updatewcs(acs_file)
    fits.setval(acs_file, keyword="NPOLFILE", value="N/A")
    updatewcs.updatewcs(acs_file)
    w1 = HSTWCS(acs_file, ext=("SCI", 1))
    w4 = HSTWCS(acs_file, ext=("SCI", 2))
    assert w1.cpdis1 is None
    assert w4.cpdis2 is None


def test_remove_d2im_distortion():
    acs_orig_file = get_filepath('j94f05bgq_flt.fits')
    current_dir = os.path.abspath(os.path.curdir)
    acs_file = get_filepath('j94f05bgq_flt.fits', current_dir)

    idctab = get_filepath('postsm4_idc.fits')
    npol_file = get_filepath('qbu16424j_npl.fits')
    d2imfile = get_filepath('new_wfc_d2i.fits ')

    try:
        os.remove(acs_file)
    except OSError:
        pass
    shutil.copyfile(acs_orig_file, acs_file)
    fits.setval(acs_file, ext=0, keyword="IDCTAB", value=idctab)
    fits.setval(acs_file, ext=0, keyword="NPOLFILE", value=npol_file)
    fits.setval(acs_file, ext=0, keyword="D2IMFILE", value=d2imfile)

    updatewcs.updatewcs(acs_file)
    fits.setval(acs_file, keyword="D2IMFILE", value="N/A")
    updatewcs.updatewcs(acs_file)
    w1 = HSTWCS(acs_file, ext=("SCI", 1))
    w4 = HSTWCS(acs_file, ext=("SCI", 2))
    assert w1.det2im1 is None
    assert w4.det2im2 is None


def test_missing_idctab():
    """ Tests that an IOError is raised if an idctab file is not found on disk."""
    acs_orig_file = get_filepath('j94f05bgq_flt.fits')
    current_dir = os.path.abspath(os.path.curdir)
    acs_file = get_filepath('j94f05bgq_flt.fits', current_dir)

    try:
        os.remove(acs_file)
    except OSError:
        pass
    shutil.copyfile(acs_orig_file, acs_file)

    fits.setval(acs_file, keyword="IDCTAB", value="my_missing_idctab.fits")
    with pytest.raises(IOError):
        updatewcs.updatewcs(acs_file)


def test_missing_npolfile():
    """ Tests that an IOError is raised if an NPOLFILE file is not found on disk."""
    acs_orig_file = get_filepath('j94f05bgq_flt.fits')
    current_dir = os.path.abspath(os.path.curdir)
    acs_file = get_filepath('j94f05bgq_flt.fits', current_dir)

    try:
        os.remove(acs_file)
    except OSError:
        pass

    shutil.copyfile(acs_orig_file, acs_file)

    fits.setval(acs_file, keyword="NPOLFILE", value="missing_npl.fits")
    with pytest.raises(IOError):
        updatewcs.updatewcs(acs_file)


def test_missing_d2imfile():
    """ Tests that an IOError is raised if a D2IMFILE file is not found on disk."""
    acs_orig_file = get_filepath('j94f05bgq_flt.fits')
    current_dir = os.path.abspath(os.path.curdir)
    acs_file = get_filepath('j94f05bgq_flt.fits', current_dir)

    try:
        os.remove(acs_file)
    except OSError:
        pass

    shutil.copyfile(acs_orig_file, acs_file)

    fits.setval(acs_file, keyword="D2IMFILE", value="missing_d2i.fits")
    with pytest.raises(IOError):
        updatewcs.updatewcs(acs_file)


def test_found_idctab():
    """ Tests the return value of apply_corrections.foundIDCTAB()."""
    acs_orig_file = get_filepath('j94f05bgq_flt.fits')
    current_dir = os.path.abspath(os.path.curdir)
    acs_file = get_filepath('j94f05bgq_flt.fits', current_dir)
    npol_file = get_filepath('qbu16424j_npl.fits')
    d2imfile = get_filepath('new_wfc_d2i.fits ')

    try:
        os.remove(acs_file)
    except OSError:
        pass

    shutil.copyfile(acs_orig_file, acs_file)
    fits.setval(acs_file, ext=0, keyword="NPOLFILE", value=npol_file)
    fits.setval(acs_file, ext=0, keyword="D2IMFILE", value=d2imfile)
    fits.setval(acs_file, keyword="IDCTAB", value="N/A")
    corrections = apply_corrections.setCorrections(acs_file)
    assert('MakeWCS' not in corrections)
    assert('TDDCor' not in corrections)
    assert('CompSIP' not in corrections)


def test_add_radesys():
    """ test that RADESYS was successfully added to headers."""
    acs_orig_file = get_filepath('j94f05bgq_flt.fits')
    current_dir = os.path.abspath(os.path.curdir)
    acs_file = get_filepath('j94f05bgq_flt.fits', current_dir)

    idctab = get_filepath('postsm4_idc.fits')
    npol_file = get_filepath('qbu16424j_npl.fits')
    d2imfile = get_filepath('new_wfc_d2i.fits ')

    try:
        os.remove(acs_file)
    except OSError:
        pass

    shutil.copyfile(acs_orig_file, acs_file)
    fits.setval(acs_file, ext=0, keyword="IDCTAB", value=idctab)
    fits.setval(acs_file, ext=0, keyword="NPOLFILE", value=npol_file)
    fits.setval(acs_file, ext=0, keyword="D2IMFILE", value=d2imfile)

    updatewcs.updatewcs(acs_file)
    for ext in [('SCI', 1), ('SCI', 2)]:
        hdr = fits.getheader(acs_file, ext)
        assert hdr['RADESYS'] == 'FK5'

def test_update_d2im_distortion():
    acs_orig_file = get_filepath('j94f05bgq_flt.fits')
    current_dir = os.path.abspath(os.path.curdir)
    acs_file = get_filepath('j94f05bgq_flt.fits', current_dir)

    idctab = get_filepath('postsm4_idc.fits')
    npol_file = get_filepath('qbu16424j_npl.fits')
    d2imfile = get_filepath('new_wfc_d2i.fits')
    newd2im = get_filepath('new_wfc_d2i.fits', current_dir)
    try:
        os.remove(acs_file)
    except OSError:
        pass
    shutil.copyfile(acs_orig_file, acs_file)
    fits.setval(acs_file, ext=0, keyword="IDCTAB", value=idctab)
    fits.setval(acs_file, ext=0, keyword="NPOLFILE", value=npol_file)
    fits.setval(acs_file, ext=0, keyword="D2IMFILE", value=d2imfile)
    updatewcs.updatewcs(acs_file)
    d2imerr1 = fits.getval(acs_file, ext=1, keyword='D2IMERR1')
    d2imerr4 = fits.getval(acs_file, ext=4, keyword='D2IMERR1')
    shutil.copyfile(d2imfile, newd2im)
    with fits.open(newd2im, mode='update') as newf:
        for ext in newf[1:]:
            ext.data = ext.data * 100

    fits.setval(acs_file, keyword="D2IMFILE", value=newd2im)
    updatewcs.updatewcs(acs_file)
    nd2imerr1 = fits.getval(acs_file, ext=1, keyword='D2IMERR1')
    nd2imerr4 = fits.getval(acs_file, ext=4, keyword='D2IMERR1')
    assert np.isclose(d2imerr1 * 100, nd2imerr1)
    assert np.isclose(d2imerr4 * 100, nd2imerr4)


def test_apply_d2im():
    from stwcs.updatewcs import apply_corrections as appc
    acs_orig_file = get_filepath('j94f05bgq_flt.fits')
    current_dir = os.path.abspath(os.path.curdir)
    fname = get_filepath('j94f05bgq_flt.fits', current_dir)
    d2imfile = get_filepath('new_wfc_d2i.fits')
    try:
        os.remove(fname)
    except OSError:
        pass
    shutil.copyfile(acs_orig_file, fname)
    fits.setval(fname, ext=0, keyword="D2IMFILE", value=d2imfile)
    fits.setval(fname, ext=0, keyword="IDCTAB", value='N/A')
    fits.setval(fname, ext=0, keyword="NPOLFILE", value='N/A')
    # If D2IMEXT does not exist, the correction should be applied
    assert appc.apply_d2im_correction(fname, d2imcorr=True)
    updatewcs.updatewcs(fname)

    # Test the case when D2IMFILE == D2IMEXT
    assert not appc.apply_d2im_correction(fname, d2imcorr=True)
    assert not appc.apply_d2im_correction(fname, d2imcorr=False)

    fits.setval(fname, ext=0, keyword='D2IMFILE', value="N/A")
    assert not appc.apply_d2im_correction(fname, d2imcorr=True)
    # No D2IMFILE keyword in primary header
    fits.delval(fname, ext=0, keyword='D2IMFILE')
    assert not appc.apply_d2im_correction(fname, d2imcorr=True)
