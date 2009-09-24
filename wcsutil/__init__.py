#from .. pywcs.sip import SIP
from pywcs import WCS, DistortionLookupTable
import pyfits
import instruments
#from .. distortion import models
from stwcs.distortion import models
import numpy as np
from pytools import fileutil
from pytools.fileutil import DEGTORAD, RADTODEG

#from .. mappings import inst_mappings, ins_spec_kw, DEGTORAD, RADTODEG, basic_wcs
import mappings
from mappings import inst_mappings, ins_spec_kw
from mappings import basic_wcs, prim_hdr_kw


__docformat__ = 'restructuredtext'
__version__ = '0.3'

class HSTWCS(WCS):
    """
    Purpose
    =======
    Create a WCS object based on the instrument.
    It has all basic WCS kw as attribbutes (set by pywcs).
    It also uses the primary and extension header to define
    instrument specific attributes.
    """
    def __init__(self, fobj=None, ext=None, instrument=None, detector=None, minerr=0.0):
        """
        :Parameters:
        `fobj`: string or PyFITS HDUList object or None
                a file name, e.g j9irw4b1q_flt.fits
                a fully qualified filename[EXTNAME,EXTNUM], e.g. j9irw4b1q_flt.fits[sci,1]
                a pyfits file object, e.g pyfits.open('j9irw4b1q_flt.fits')
        `ext`:  int or None
                extension number
                if ext==None, it is assumed the data is in the primary hdu
        `instrument`: string
                one of 'ACS', 'NICMOS', 'WFPC2', 'STIS', 'WFC3'
                Used only to define a default HSTWCS object, when fobj==None
        `detector`:  string
                for example 'WFC'
                If instrument and detector parameters are given, a default HSTWCS
                instrument is created. Used only with fobj==None
        `minerr`: float
                minimum value a distortion correction must have in order to be applied.
                If CPERRja, CQERRja are smaller than minerr, the corersponding
                distortion is not applied.
        """

        self.inst_kw = ins_spec_kw
        self.minerr = minerr
        if instrument == None:
            filename, hdr0, ehdr, phdu = self.parseInput(f=fobj, ext=ext)
            self.filename = filename
            self.setHDR0kw(hdr0, ehdr)
            WCS.__init__(self, ehdr, fobj=phdu, minerr=minerr)
            # If input was a pyfits HDUList object, it's the user's
            # responsibility to close it, otherwise, it's closed here.
            if not isinstance(fobj, pyfits.HDUList):
                phdu.close()

            self.instrument ='DEFAULT'
            self.setInstrSpecKw(hdr0, ehdr)
            self.setPscale()
            self.setOrient()
            self.readIDCCoeffs(ehdr)
            extname = ehdr.get('EXTNAME', "")
            extnum = ehdr.get('EXTVER', None)
            self.extname = (extname, extnum)
        else:
            # create a default HSTWCS object based on
            #instrument and detector parameters
            self.instrument = instrument
            self.detector = detector
            self.setInstrSpecKw()

    def parseInput(self, f=None, ext=None):
        if isinstance(f, str):
            # create an HSTWCS object from a filename

            if ext != None:
                filename = f
                if isinstance(ext,tuple):
                    if ext[0] == '':
                        extnum = ext[1] # handle ext=('',1)
                    else:
                        extnum = ext
                else:
                    extnum = int(ext)
            elif ext == None:
                filename, ext = fileutil.parseFilename(f)
                ext = fileutil.parseExtn(ext)
                if ext[0] == '':
                    extnum = int(ext[1]) #handle ext=('',extnum)
                else:
                    extnum = ext
            phdu = pyfits.open(filename)
            hdr0 = pyfits.getheader(filename)
            try:
                ehdr = pyfits.getheader(filename, ext=extnum)
            except (IndexError,KeyError):
                print 'Unable to get extension.', extnum
                raise

        elif isinstance(f, pyfits.HDUList):
            phdu = f
            if ext == None:
                extnum = 0
            else:
                extnum = ext
            ehdr = f[extnum].header
            hdr0 = f[0].header
            filename = hdr0.get('FILENAME', "")

        return filename, hdr0, ehdr, phdu

    def setHDR0kw(self, primhdr, ehdr):
        # Set attributes from kw defined in the primary header.
        self.instrument = primhdr.get('INSTRUME', None)
        self.offtab = primhdr.get('OFFTAB', None)
        self.idctab = primhdr.get('IDCTAB', None)
        self.date_obs = primhdr.get('DATE-OBS', None)
        self.ra_targ = primhdr.get('RA_TARG', None)
        self.dec_targ = primhdr.get('DEC_TARG', None)
        self.det2imfile = primhdr.get('D2IMFILE', None)
        self.det2imext = ehdr.get('D2IMEXT', None)
        try:
            self.pav3 = primhdr['PA_V3']

        except KeyError:
            print 'Keyword PA_V3 not found in primary header.'
            print 'This is typical for some old files. Please retrieve the files from the archive again.'
            #print 'Quitting ...'
            #raise



    def readIDCCoeffs(self, header):
        """
        Reads in first order IDCTAB coefficients if present in the header
        """
        coeffs = ['ocx10', 'ocx11', 'ocy10', 'ocy11', 'idcscale']
        for c in coeffs:
            self.__setattr__(c, header.get(c, None))

    def setInstrSpecKw(self, prim_hdr=None, ext_hdr=None):
        # Based on the instrument kw creates an instalnce of an instrument WCS class
        # and sets attributes from instrument specific kw
        if self.instrument in inst_mappings.keys():
            inst_kl = inst_mappings[self.instrument]
            inst_kl = instruments.__dict__[inst_kl]
            insobj = inst_kl(prim_hdr, ext_hdr)
            for key in self.inst_kw:
                try:
                    self.__setattr__(key, insobj.__getattribute__(key))
                except AttributeError:
                    # Some of the instrument's attributes are recorded in the primary header and
                    # were already set, (e.g. 'DETECTOR'), the code below is a check for that case.
                    if not self.__getattribute__(key):
                        raise
                    else:
                        pass
        else:
            raise KeyError, "Unsupported instrument - %s" %self.instrument

    def setPscale(self):
        # Calculates the plate scale from the cd matrix

        cd11 = self.wcs.cd[0][0]
        cd21 = self.wcs.cd[1][0]
        self.pscale = np.sqrt(np.power(cd11,2)+np.power(cd21,2)) * 3600.

    def setOrient(self):
        # Recompute ORIENTAT
        cd12 = self.wcs.cd[0][1]
        cd22 = self.wcs.cd[1][1]
        self.orientat = RADTODEG(np.arctan2(cd12,cd22))

    def updatePscale(self, pscale):
        """Given a plate scale, update the CD matrix"""
        old_pscale = self.pscale
        self.pscale = pscale
        self.wcs.cd = self.wcs.cd * pscale/old_pscale
        self.naxis1 = self.naxis1 * old_pscale/ pscale
        self.naxis2 = self.naxis2 * old_pscale/ pscale
        self.wcs.crpix = self.wcs.crpix *old_pscale/pscale

    def updateOrient(self, orient):
        """Given n angle update the CD matrix"""
        if self.orientat == orient:
            return
        old_orient = self.orientat
        self.orientat = orient
        angle = fileutil.DEGTORAD(orient)
        cd11 = -np.cos(angle)
        cd12 = np.sin(angle)
        cd21 = cd12
        cd22 = -cd11
        cdmat = np.array([[cd11, cd12],[cd21,cd22]])
        self.wcs.cd = cdmat * self.pscale/3600


    def readModel(self, update=False, header=None):
        """
        Purpose
        =======
        Read distortion model from idc table.
        Save some of the information as kw needed for interpreting the distortion
        If header is provided and update is True, some IDC model kw
        will be recorded in the header.
        """
        if self.idctab == None or self.date_obs == None:
            print 'idctab or date_obs not available\n'
            self.idcmodel = None
            return
        if self.filter1 ==None and self.filter2 == None:
            'No filter information available\n'
            self.idcmodel = None
            return

        self.idcmodel = models.IDCModel(self.idctab,
                    chip=self.chip, direction='forward', date=self.date_obs,
                    filter1=self.filter1, filter2=self.filter2,
                    offtab=self.offtab, binned=self.binned)

        if update:
            if header==None:
                print 'Update header with IDC model kw requested but header was not provided\n.'
                return
            else:
                self.updatehdr(header)
    
    
    def restore(self, header=None):
        """
        Restore a WCS archive in memory and update the WCS object.
        Restored are the basic WCS kw as well as pscale and orient.
        """
        from pywcs import Wcsprm
        backup = {}
        if header == None:
            print 'Need a valid header in order to restore archive\n'
            return

        for k in basic_wcs:
            try:
                nkw = ('O'+k)[:7]
                backup[k] = header[nkw]
            except KeyError:
                pass
        if backup == {}:
            print 'No archive was found\n'
            return
        cdl=pyfits.CardList()
        for item in backup.items():
            card = pyfits.Card(key=item[0], value=item[1])
            cdl.append(card)

        h = pyfits.Header(cdl)
        wprm = Wcsprm("".join([str(x) for x in h.ascardlist()]))
        self.wcs = wprm
        self.setPscale()
        self.setOrient()

    def updatehdr(self, ext_hdr, newkeywords=None):
        #kw2add : OCX10, OCX11, OCY10, OCY11
        # record the model in the header for use by pydrizzle
        ext_hdr.update('OCX10', self.idcmodel.cx[1,0])
        ext_hdr.update('OCX11', self.idcmodel.cx[1,1])
        ext_hdr.update('OCY10', self.idcmodel.cy[1,0])
        ext_hdr.update('OCY11', self.idcmodel.cy[1,1])
        ext_hdr.update('IDCSCALE', self.idcmodel.refpix['PSCALE'])
        ext_hdr.update('IDCTHETA', self.idcmodel.refpix['THETA'])
        ext_hdr.update('IDCXREF', self.idcmodel.refpix['XREF'])
        ext_hdr.update('IDCYREF', self.idcmodel.refpix['YREF'])
        ext_hdr.update('IDCV2REF', self.idcmodel.refpix['V2REF'])
        ext_hdr.update('IDCV3REF', self.idcmodel.refpix['V3REF'])

    def printwcs(self):
        print 'WCS Keywords\n'
        print 'CD_11  CD_12: %r %r' % (self.wcs.cd[0,0],  self.wcs.cd[0,1])
        print 'CD_21  CD_22: %r %r' % (self.wcs.cd[1,0],  self.wcs.cd[1,1])
        print 'CRVAL    : %r %r' % (self.wcs.crval[0], self.wcs.crval[1])
        print 'CRPIX    : %r %r' % (self.wcs.crpix[0], self.wcs.crpix[1])
        print 'NAXIS    : %d %d' % (self.naxis1, self.naxis2)
        print 'Plate Scale : %r' % self.pscale
        print 'ORIENTAT : %r' % self.orientat

def help():
    print 'How to create an HSTWCS object:\n\n'
    print """ \
    1. Using a pyfits HDUList object and an extension number \n
    Example:\n

    fobj = pyfits.open('some_file.fits')\n
    w = wcsutil.HSTWCS(fobj, 3)\n\n

    2. Create an HSTWCS object using a qualified file name. \n
    Example:\n
    w = wcsutil.HSTWCS('j9irw4b1q_flt.fits[sci,1]')\n\n

    3. Create an HSTWCS object using a file name and an extension number. \n
    Example:\n
    w = wcsutil.HSTWCS('j9irw4b1q_flt.fits', ext=2)\n\n

    4. Create a template HSTWCS object for an instrument/detector combination.]n
    Example:\n
    w = wcsutil.HSTWCS(instrument='ACS', detector='WFC'\n\n
    """

