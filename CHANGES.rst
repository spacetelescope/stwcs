Unreleased
----------
- When the GSC web service returns no RA and DEC offsets from GAIA (0.), use default FITS WCS values.

1.7.4 (2025-01-07)
------------------

- Manual scalar promotion of refpix['PSCALE'] from float32
  to float64 to avoid future Numpy 2.0 issues [#206].

- Pin numpy min version greater than 2.0 [#214] 


1.7.3 (2024-05-06)
------------------
  
- Avoid Exception for some new data. [#189]
  
- Convert WCSNAME to hash for HLET filename. [#183]
  
- Pin astropy min version to 5.0.4 [#191]
  

1.7.2 (2021-11-29)
-----------------

Bug Fixes
^^^^^^^^^

- Fixed bug affecting data with no astrometry database solutions so
  that they can now be processed successfully. [#189]

1.7.1 (2021-10-20)
------------------
- Filenames for new headerlet files based on updated/new
  IDCTABs now based on hash of WCSNAME instead of WCSNAME itself. [#183]


1.7.0 (2021-08-10)
------------------

- Make the code more robust when defining ``extname`` attribute of a WCS for
  multi-extension FITS files when their primary header contains ``'EXTVER'``
  keyword but not the ``'EXTNAME'`` keywords. [#179]

- Fix a bug in ``stwcs.wcsutil.altwcs.archive_wcs`` due to which
  a ``str`` type ``ext`` would raise an exception. [#177]

- Modify ``updatewcs`` to insure that the distortion keywords remain
  in the same part of the primary header as the rest of the calibration
  reference file names. [#175]

- Update ``build_reference_wcs`` to accept list of filenames as well
  single filenames or a single ``HDUList`` object as input. [#173]

- Add support for computing new `a priori` WCS solutions when running
  ``updatewcs``. [#170]

- Add new parameter to ``updatewcs`` to limit the astrometry database
  WCSs appended to the file to only those based on the same IDCTAB as
  specified in the image header. [#170]

- Add new parameter to ``updatewcs`` to specify whether or not to
  remove duplicate headerlet extensions from the FITS file when updating
  the image using the astrometry database. [#172]


1.6.1 (2020-12-09)
------------------

- Fix a bug in ``wcsutil.wcscorr.init_wcscorr()`` that would result in crash
  when run on drizzled images that do not have OPUS WCS (WCS key 'O'). [#165]

- Make WCS and headerlet name comparisons case-insensitive when applying
  headerlets. [#163]

- Deprecate ``accuracy`` argument in ``all_world2pix`` and replace it with
  ``tolerance`` in order to have compatible function signature with
  ``astropy``'s ``all_world2pix()``'. [#166]


1.6.0 (2020-07-16)
------------------

- Prevent same WCS from being updated from astrometry database headerlet. [#111]

- Revise how headerlets are applied as primary WCS [#122]

- Add comments to WCS keywords for ``IDC*``, ``OCX``, ``OCY``, and alternate
  WCS. [#127]

- Primary WCS keywords no longer contain trailing spaces which, in the past,
  lead to duplicate keywords in image headers. [#131]

- AltWCS: Non-standard WCS keywords defined by, e.g., ``TweakReg`` are now
  archived, read, and restored together with standard WCS keywords. [#135]

- WCS keywords that are not specific to an alternate WCS (such as ``MJDREF``)
  are no longer archived. [#141]

- Fix a bug in ``_restore()`` due to which ``TDDALPHA`` and ``TDDBETA`` were
  not reset to ``0.0`` when restoring wrom the ``'OPUS'`` WCS with key ``'O'``. [#141]

- Increased reliability of detecting alternate WCS in image headers. [#141]

- Abbreviated comments for ``D2IM``, ``CPDIS`` lookup table distortions. [#142]

- Add comments to ``WCSNAMEa`` header keywords. [#143]

- Do not archive primary WCS when applying headerlet as primary when existing
  primary WCS was created from the same headerlet that is being applied, based
  on the value of ``'HDRNAME'`` keyword. [#147]

- Added ``archive_wcs()`` function with more advanced control over how to
  archive primary WCS and with more capabilities, such as auto-renaming. [#153]

- Remove HST-specific WCS keywords from SCI extensions that, following HST
  convention, are supposed to be present in primary header only. [#159]

- Do not archive primary WCS when applying headerlet as primary if distortion
  model changes. [#161]


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
