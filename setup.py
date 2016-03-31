#!/usr/bin/env python
import relic.release
from glob import glob
from numpy import get_include as np_include
from setuptools import setup, find_packages, Extension


version = relic.release.get_info()
relic.release.write_template(version, 'lib/stwcs')

setup(
    name = 'stwcs',
    version = version.pep386,
    author = 'Nadia Dencheva, Warren Hack',
    author_email = 'help@stsci.edu',
    description = 'Recomputes the WCS of an HST observation and puts all distortion',
    url = 'https://github.com/spacetelescope/stwcs',
    classifiers = [
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    install_requires = [
        'astropy',
        'nose',
        'numpy',
        'stsci.tools',
    ],
    package_dir = {
        '': 'lib',
    },
    packages = find_packages('lib'),
    package_data = {
        'stwcs/gui': ['*.help'],
        'stwcs/gui/pars': ['*'],
        'stwcs/gui/htmlhelp': ['*'],
    },
)
