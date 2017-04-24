1.3 (Unreleased)
----------------

- The ``clobber`` parameter in `Headerlet.tofile()`` was replaced with
  ``overwrite``. [#24]
- Fixed a python compatibility bug. [#30]
-

1.2.5 (2016-12-20)
----------------

- updatewcs() now reads all extension immediately after opening a file
  to fix a problem after astropy implemented fits lazy loading. [#21]

- Fixed a bug in updating the D2IM correction in a science file when the
  a new distortion file was supplied through D2IMFILE keyword. [#22]

- Warning messages from astropy.wcs are filtered out when they are not relevant. [#31]

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
