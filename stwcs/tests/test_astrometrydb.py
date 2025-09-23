import shutil
import os

import pytest
import requests

from numpy.testing import assert_allclose

from astropy.io import fits
from astropy.io.fits import diff
from .. import updatewcs
from ..updatewcs import astrometry_utils

from . import data
data_path = os.path.split(os.path.abspath(data.__file__))[0]

# set environment variable to insure Exceptions are raised
os.environ['RAISE_PIPELINE_ERRORS'] = 'True'


def get_filepath(filename, directory=data_path):
    return os.path.join(directory, filename)


class TestAstrometryDB:

    def setup_class(self):
        os.environ['ASTROMETRY_STEP_CONTROL'] = 'On'
        self.obsname = 'j94f05bgq_flt.fits'
        refname = self.obsname.replace('_flt', '_flt_rdb')
        acs_orig_file = get_filepath(self.obsname)
        current_dir = os.path.abspath(os.path.curdir)
        self.ref_file = get_filepath(refname, current_dir)
        self.acs_file = get_filepath(self.obsname, current_dir)

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

        shutil.copyfile(self.acs_file, self.ref_file)

        updatewcs.updatewcs(self.acs_file, use_db=False)

    @pytest.mark.skip("Need to understand why this fails and how it's supposed to work.")
    def test_default(self):
        """
        Sanity check: Insure it will run at all in default mode
        """
        updatewcs.updatewcs(self.ref_file)
        adb = astrometry_utils.AstrometryDB()
        adb.updateObs(self.acs_file)
        # at this point self.acs_file == self.ref_file if all worked...
        acs = fits.open(self.acs_file)
        ref = fits.open(self.ref_file)
        report = diff.HDUDiff(acs[1], ref[1], ignore_keywords=['HDRNAME', 'HDRNAMEB']).report()
        assert "No differences found" in report

    def test_new_obs(self, caplog):
        """
        A simple sanity check that first time processing will not crash
        since that observation name will not be found in the database.
        """
        new_obsname = self.obsname.replace('j94', 'a94')
        shutil.copyfile(self.acs_file, new_obsname)
        fits.setval(new_obsname, ext=0, keyword="rootname", value="a94f05bgq")
        adb = astrometry_utils.AstrometryDB(perform_step=True)
        adb.updateObs(new_obsname)
        assert adb.new_observation
        os.remove(new_obsname)  # remove intermediate test file
        del adb

    def test_no_offsets(self):
        """ HLA-1541"""
        new_obsname = 'j8di67a2q_flt.fits'
        shutil.copyfile(self.acs_file, new_obsname)
        fits.setval(new_obsname, keyword="rootname", value="j8di67a2q")
        offsets = astrometry_utils.find_gsc_offset('j8di67a2q_flt.fits')
        expected = {
            'delta_x': 0.0,
            'delta_y': 0.0,
            'delta_ra': 0.0,
            'delta_dec': 0.0,
            'dGSinputDEC': 0.0,
            'dGSinputRA': 0.0,
            'dGSoutputDEC': 0.0,
            'dGSoutputRA': 0.0,
            'roll': 0.0,
            'scale': 1.0,
            'catalog': None,
            }
        # Do not compare the WCS
        offsets.pop("expwcs")
        offsets.pop("message")
        assert expected == offsets
        os.remove(new_obsname)

        # mock a request response of False, code 403
        os.environ['GSSS_WEBSERVICES_URL']="http://test.com"
        serviceUrl='http://test.com/GSCConvert/GSCconvert.aspx?IPPPSSOOT=a8di67a2q'
        rawcat= requests.get(serviceUrl)
        assert not rawcat.ok
        offsets = astrometry_utils.find_gsc_offset(self.acs_file)
        offsets.pop("expwcs")
        offsets.pop("message")
        assert expected == offsets
        del os.environ['GSSS_WEBSERVICES_URL']

    def test_success_offsets(self):
        """Test successfully retrieving offsets from the GSC service."""
        expected = {'delta_ra': 0.00163279,
                    'delta_dec': -0.00024635,
                    'roll': 0.05507061,
                    'scale': 0.99852027,
                    'dGSinputRA': 4.870375,
                    'dGSinputDEC': -72.181833,
                    'dGSoutputRA': 4.87200779,
                    'dGSoutputDEC': -72.18207935,
                    'catalog': 'GSC240',
                    'delta_x': -26.195118548980645,
                    'delta_y': -30.37243670094358
                   }

        offsets = astrometry_utils.find_gsc_offset(self.acs_file)
        offsets.pop("expwcs")
        offsets.pop("message")
        assert expected.pop('catalog') == offsets.pop('catalog')
        for item in list(expected.keys())[::-1]:
            # can't compare directly because computed results are slightly different on
            # different OSs.
            assert_allclose(expected.pop(item), offsets.pop(item))


def test_db_connection():
    adb = astrometry_utils.AstrometryDB()
    assert adb.available
    del adb


def test_db_raise_errors_user_override():
    # should not raise an exception since raise_errors is False
    db = astrometry_utils.AstrometryDB(url="bad_link/", raise_errors=False)
    assert db.available is False
    del db


@pytest.mark.xfail
def test_db_raise_errors_true():
    # tests that environment variable RAISE_PIPELINE_ERRORS raises an exception
    db = astrometry_utils.AstrometryDB(url="bad_link/")
    del db


@pytest.mark.xfail
def test_db_timeout():
    # Testing mode uses short timeout to simulate timeout. raise_errors should
    # be set to True with the environment variable.
    db = astrometry_utils.AstrometryDB()
    db.isAvailable(testing=True)
    del db
