.. Headerlet Definition TSR master file

TSR 2012-02: Definition of a Headerlet and Its Role in Updating WCS Information
===============================================================================

   Abstract::
   Authors: Warren Hack, Nadezhda Dencheva
   Date: 22 Oct 2012

   A headerlet is a self-consistent representation of a single WCS solution for a single
   exposure complete with all distortion information. FITS is the data
   storage format currently supported. It has no observational data
   which makes it relatively small and light to distribute.
   It is, essentially, a mechanism for encapsulating WCS information
   which can later be used to update the WCS of a science file and allows
   improved astrometric solutions to be stored and passed around easily.
   The HST archive is expected to start accepting headerlets for HST data soon.
   However the implementaiton is not HST specific and FITS WCS standard, as well as
   all WCS conventions implemented in pywcs are supported.
   This report describes the format and contents of a headerlet
   along with the software implementation and methods for creating headerlets and using them
   to update the WCS of a science observation.

Contents:

.. toctree::
   :maxdepth: 2

   hdrlet_tsr
