#!/usr/bin/env python
import os
import subprocess
import sys
from glob import glob
from numpy import get_include as np_include
from setuptools import setup, find_packages, Extension


if os.path.exists('relic'):
    sys.path.insert(1, 'relic')
    import relic.release
else:
    try:
        import relic.release
    except ImportError:
        try:
            subprocess.check_call(['git', 'clone',
                'https://github.com/jhunkeler/relic.git'])
            sys.path.insert(1, 'relic')
            import relic.release
        except subprocess.CalledProcessError as e:
            print(e)
            exit(1)


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
