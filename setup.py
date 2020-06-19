#!/usr/bin/env python
import sys
from setuptools import setup, find_packages
from configparser import ConfigParser

if sys.version_info < (3, 6):
    error = """
    STWCS supports Python 3.6 and above.
    """
    sys.exit(error)

conf = ConfigParser()
conf.read(['setup.cfg'])

# Get some config values
metadata = dict(conf.items('metadata'))
PACKAGENAME = metadata.get('package_name', 'stwcs')
DESCRIPTION = metadata.get('description', '')
AUTHOR = metadata.get('author', 'STScI')
AUTHOR_EMAIL = metadata.get('author_email', 'help@stsci.edu')

DOCS_REQUIRE = ["sphinx",
                "sphinx-automodapi",
                "sphinx-rtd-theme",
                'sphinx-automodapi',
                ]

TESTS_REQUIRE = ["pytest"]

setup(
    name = PACKAGENAME,
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
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    python_requires='>=3.6',
    install_requires = [
        'astropy',
        'numpy',
        'stsci.tools>=3.6',
        'requests',
        'lxml'
    ],
    packages = find_packages(),
    extras_require={
        'docs': DOCS_REQUIRE,
        'test': TESTS_REQUIRE,
        },
    tests_require = TESTS_REQUIRE,
    package_data = {
        'stwcs/gui': ['*.help'],
        'stwcs/gui/pars': ['*'],
        'stwcs/gui/htmlhelp': ['*'],
    },
)
