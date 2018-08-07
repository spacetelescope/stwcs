#!/usr/bin/env python
import os
import pkgutil
import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
from subprocess import check_call, CalledProcessError


if sys.version_info < (3, 5):
    error = """
    STWCS 1.4+ does not support Python 2.x, 3.0, 3.1, 3.2, 3.3 or 3.4.
    Beginning with STWCS 1.4.0, Python 3.5 and above is required.

    This may be due to an out of date pip.

    Make sure you have pip >= 9.0.1.

    """
    sys.exit(error)


if not pkgutil.find_loader('relic'):
    relic_local = os.path.exists('relic')
    relic_submodule = (relic_local and
                       os.path.exists('.gitmodules') and
                       not os.listdir('relic'))
    try:
        if relic_submodule:
            check_call(['git', 'submodule', 'update', '--init', '--recursive'])
        elif not relic_local:
            check_call(['git', 'clone', 'https://github.com/spacetelescope/relic.git'])

        sys.path.insert(1, 'relic')
    except CalledProcessError as e:
        print(e)
        exit(1)

import relic.release

version = relic.release.get_info()
relic.release.write_template(version, 'stwcs')

try:
    from distutils.config import ConfigParser
except ImportError:
    from configparser import ConfigParser

conf = ConfigParser()
conf.read(['setup.cfg'])

# Get some config values
metadata = dict(conf.items('metadata'))
PACKAGENAME = metadata.get('package_name', 'stwcs')
DESCRIPTION = metadata.get('description', '')
AUTHOR = metadata.get('author', 'STScI')
AUTHOR_EMAIL = metadata.get('author_email', 'help@stsci.edu')

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['stwcs/tests']
        self.test_suite = True

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

setup(
    name = PACKAGENAME,
    version = version.pep386,
    author = AUTHOR,
    author_email = AUTHOR_EMAIL,
    description = DESCRIPTION,
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
        'numpy',
        'stsci.tools',
        'requests',
        'lxml'
    ],
    packages = find_packages(),
    tests_require = ['pytest'],
    package_data = {
        'stwcs/gui': ['*.help'],
        'stwcs/gui/pars': ['*'],
        'stwcs/gui/htmlhelp': ['*'],
    },
    cmdclass = {"test": PyTest}
)
