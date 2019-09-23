1.5.3 (2019-09-23)
------------------

- Add ``stdout`` to logging handlers. [#108]

- Correct the logic for replacing headerlets with one which has a different 
  distortion model. [#109]


1.5.2 (2019-08-26)
------------------

- Correct a problem where best WCS solutions were not applied
  with repeated reruns of ``updatewcs``. [#95]

1.5.1 (2019-07-03)
------------------

- Improve Travis test matrices to include testing with public and dev
  versions of dependencies. [#87]

- Fixed a bug in converting a ``PC`` to a ``CD`` matrix. [#77]

- Filter out expected warnings from tests. [#78]

- Allow ``updatewcs`` to be called with ``HDUList`` objects as input. [#80]

- Update the XML parser for the astrometry database and switch the default to use
  the MAST TEST server which is publicly accessible. [#74]

- Gracefully ignore when the astrometry database returns an empty result for
  an image. [#84]

- Insure updatewcs works with WFPC2 data while supporting HDUList inputs. [#85]

1.4.2(2018-08-07)
-----------------

- Fix a problem with the build script. [#64]

1.4.1(2018-08-06)
-----------------
- Fix a bug in restoring headerlets from a HeaderletHDU to a SCI extension [#60]

- Fix logic to gracefully handle lack of any WCS solutions from the
  astrometry database for an exposure.  [#62]

1.4.0(2018-01-22)
-----------------

- Fix a bug in creating headerlets from a I/O stream. [#39]

- Added an interface to a new astrometry database which will
  contain new astrometrically-accurate solutions for the pointing
  for HST observations. [#40]

1.3.2 (20017-07-05)
-------------------

- The ``clobber`` parameter in `Headerlet.tofile()`` was replaced with
  ``overwrite``. [#24]

- Fixed a python compatibility bug. [#30]

- Warning messages from astropy.wcs are filtered out when they are not relevant. [#31]


1.2.5 (2016-12-20)
----------------

- updatewcs() now reads all extension immediately after opening a file
  to fix a problem after astropy implemented fits lazy loading. [#21]

- Fixed a bug in updating the D2IM correction in a science file when the
  a new distortion file was supplied through D2IMFILE keyword. [#22]

1.2.4 (2016-10-27)
------------------

- Fix a bug in removing LookupTable distortion. [#3]

- Fix a ``KeyError`` crash in applying headerlets. [#9]

- Catch ``KeyError`` when deleting header keywords. [#13]

- Fix for ``REFFRAME = OTHER``. [#14]

- In cases when the warning is expected, catch INFO messages
  coming from `astropy.wcs`. [#16]


1.2.3 (2016-07-13)
------------------

- Move to github.
