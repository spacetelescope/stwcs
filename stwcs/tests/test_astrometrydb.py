import shutil
import os

import pytest
<<<<<<< HEAD
=======

>>>>>>> acb99e7 (processing with no offsets from the GSC)
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

<<<<<<< HEAD
    @pytest.mark.skip("Need to understand why this fails and how it's supposed to work.")
=======
    def test_db_connection(self):

        adb = astrometry_utils.AstrometryDB()
        assert adb.available
        del adb

<<<<<<< HEAD
>>>>>>> acb99e7 (processing with no offsets from the GSC)
=======
    @pytest.mark.skip("Need to understand why this fails and how it's supposed to work.")
>>>>>>> 907a89c (add and skip another test)
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
<<<<<<< HEAD
        assert "No differences found" in report
=======
        print(report)
        assert report == ""
>>>>>>> 907a89c (add and skip another test)

    @pytest.mark.skip("Need to understand why this fails and how it's supposed to work.")
    def test_new_obs(self):
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

<<<<<<< HEAD

def test_db_connection():

    adb = astrometry_utils.AstrometryDB()
    adb.isAvailable()
    assert adb.available
    del adb
=======
    def test_no_offset(self):
        """ HLA-1541"""
        new_obsname = 'j8di67a2q_flt.fits'
        shutil.copyfile(self.acs_file, new_obsname)
        fits.setval(new_obsname, keyword="rootname", value="j8di67a2q")
        offsets = astrometry_utils.find_gsc_offset('j8di67a2q_flt.fits')
        expected = {
            'delta_x': 0.0,
            'delta_y': 0.0,
            'roll': 0.0,
            'scale': 1.0,
            'delta_ra': 0.0,
            'delta_dec': 0.0,
            'catalog': None
            }
        # Do not compare the WCS
        offsets.pop("expwcs")
        assert expected == offsets
        os.remove(new_obsname)
>>>>>>> acb99e7 (processing with no offsets from the GSC)
