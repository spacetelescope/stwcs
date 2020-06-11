import shutil
import os

from astropy.io import fits
from .. import updatewcs
from ..updatewcs import astrometry_utils

from . import data
data_path = os.path.split(os.path.abspath(data.__file__))[0]

# set environment variable to insure Exceptions are raised
os.environ['RAISE_PIPELINE_ERRORS'] = 'True'


def get_filepath(filename, directory=data_path):
    return os.path.join(directory, filename)


class TestAstrometryDB(object):

    def setup_class(self):
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

    def test_db_connection(self):

        adb = astrometry_utils.AstrometryDB()
        adb.isAvailable()
        del adb

    def test_default(self):
        """
        Sanity check: Insure it will run at all in default mode
        """

        updatewcs.updatewcs(self.ref_file)

        adb = astrometry_utils.AstrometryDB()
        adb.updateObs(self.acs_file)
        # at this point self.acs_file == self.ref_file if all worked...

    def test_new_obs(self):
        """
        A simple sanity check that first time processing will not crash
        since that observation name will not be found in the database.
        """
        new_obsname = self.obsname.replace('j94', 'a94')
        shutil.copyfile(self.acs_file, new_obsname)
        adb = astrometry_utils.AstrometryDB()
        adb.updateObs(new_obsname)

        os.remove(new_obsname)  # remove intermediate test file
        del adb
