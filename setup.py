from distutils.core import setup
from os.path import join
# PyFITS
try:
    import pyfits
except ImportError:
    print "WARNING: PyFITS must be installed to use hstwcs."
    print "         Since this is not a build-time dependency, the build will proceed."

# PyWCS
try:
    import pywcs
except ImportError:
    print "WARNING: PyWCS must be installed to use hstwcs."
    print "         Since this is not a build-time dependency, the build will proceed."

setup(name="hstwcs",
      version="0.1",
      description="HST WCS Corrections",
      packages=['hstwcs', 'hstwcs/updatewcs', 'hstwcs/wcsutil', 'hstwcs/distortion'],
      package_dir={'hstwcs':'lib', 'hstwcs/updatewcs': 'updatewcs',
                    'hstwcs/wcsutil': 'wcsutil', 'hstwcs/distortion': 'distortion'}
      )
