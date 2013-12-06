from __future__ import division # confidence high

import os.path
from pywcs import WCS
import pyfits
import instruments
from stwcs.distortion import models, coeff_converter
import altwcs
import numpy as np
from stsci.tools import fileutil

import getinput
import mappings
from mappings import inst_mappings, ins_spec_kw
from mappings import basic_wcs


__docformat__ = 'restructuredtext'

#
#### Utility functions copied from 'updatewcs.utils' to avoid circular imports
#
def extract_rootname(kwvalue,suffix=""):
    """ Returns the rootname from a full reference filename

        If a non-valid value (any of ['','N/A','NONE','INDEF',None]) is input,
            simply return a string value of 'NONE'

        This function will also replace any 'suffix' specified with a blank.
    """
    # check to see whether a valid kwvalue has been provided as input
    if kwvalue.strip() in ['','N/A','NONE','INDEF',None]:
        return 'NONE' # no valid value, so return 'NONE'

    # for a valid kwvalue, parse out the rootname
    # strip off any environment variable from input filename, if any are given
    if '$' in kwvalue:
        fullval = kwvalue[kwvalue.find('$')+1:]
    else:
        fullval = kwvalue
    # Extract filename without path from kwvalue
    fname = os.path.basename(fullval).strip()

    # Now, rip out just the rootname from the full filename
    rootname = fileutil.buildNewRootname(fname)

    # Now, remove any known suffix from rootname
    rootname = rootname.replace(suffix,'')
    return rootname.strip()

def build_default_wcsname(idctab):

    idcname = extract_rootname(idctab,suffix='_idc')
    wcsname = 'IDC_' + idcname
    return wcsname

#
#### HSTWCS Class definition
#
class HSTWCS(WCS):

    def __init__(self, fobj=None, ext=None, minerr=0.0, wcskey=" "):
        """
        Create a WCS object based on the instrument.

        In addition to basic WCS keywords this class provides
        instrument specific information needed in distortion computation.

        Parameters
        ----------
        fobj: string or PyFITS HDUList object or None
                a file name, e.g j9irw4b1q_flt.fits
                a fully qualified filename[EXTNAME,EXTNUM], e.g. j9irw4b1q_flt.fits[sci,1]
                a pyfits file object, e.g pyfits.open('j9irw4b1q_flt.fits'), in which case the
                user is responsible for closing the file object.
        ext:  int, tuple or None
                extension number
                if ext is tuple, it must be ("EXTNAME", EXTNUM), e.g. ("SCI", 2)
                if ext is None, it is assumed the data is in the primary hdu
        minerr: float
                minimum value a distortion correction must have in order to be applied.
                If CPERRja, CQERRja are smaller than minerr, the corersponding
                distortion is not applied.
        wcskey: str
                A one character A-Z or " " used to retrieve and define an
                alternate WCS description.
        """

        self.inst_kw = ins_spec_kw
        self.minerr = minerr
        self.wcskey = wcskey

        if fobj is not None:
            filename, hdr0, ehdr, phdu = getinput.parseSingleInput(f=fobj,
                                                                   ext=ext)
            self.filename = filename
            instrument_name = hdr0.get('INSTRUME', 'DEFAULT')
            if instrument_name in ['IRAF/ARTDATA','',' ','N/A']:
                self.instrument = 'DEFAULT'
            else:
                self.instrument = instrument_name
            WCS.__init__(self, ehdr, fobj=phdu, minerr=self.minerr,
                         key=self.wcskey)
            # If input was a pyfits HDUList object, it's the user's
            # responsibility to close it, otherwise, it's closed here.
            if not isinstance(fobj, pyfits.HDUList):
                phdu.close()
            self.setInstrSpecKw(hdr0, ehdr)
            self.readIDCCoeffs(ehdr)
            extname = ehdr.get('EXTNAME', '')
            extnum = ehdr.get('EXTVER', None)
            self.extname = (extname, extnum)
        else:
            # create a default HSTWCS object
            self.instrument = 'DEFAULT'
            WCS.__init__(self, minerr=self.minerr, key=self.wcskey)
            self.pc2cd()
            self.setInstrSpecKw()
        self.setPscale()
        self.setOrient()

    def readIDCCoeffs(self, header):
        """
        Reads in first order IDCTAB coefficients if present in the header
        """
        coeffs = ['ocx10', 'ocx11', 'ocy10', 'ocy11', 'idcscale']
        for c in coeffs:
            self.__setattr__(c, header.get(c, None))

    def setInstrSpecKw(self, prim_hdr=None, ext_hdr=None):
        """
        Populate the instrument specific attributes:

        These can be in different headers but each instrument class has knowledge
        of where to look for them.

        Parameters
        ----------
        prim_hdr: pyfits.Header
                  primary header
        ext_hdr:  pyfits.Header
                  extension header

        """
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
        """
        Calculates the plate scale from the CD matrix
        """
        try:
            cd11 = self.wcs.cd[0][0]
            cd21 = self.wcs.cd[1][0]
            self.pscale = np.sqrt(np.power(cd11,2)+np.power(cd21,2)) * 3600.
        except AttributeError:
            if self.wcs.has_cd():
                print "This file has a PC matrix. You may want to convert it \n \
                to a CD matrix, if reasonable, by running pc2.cd() method.\n \
                The plate scale can be set then by calling setPscale() method.\n"
            self.pscale = None

    def setOrient(self):
        """
        Computes ORIENTAT from the CD matrix
        """
        try:
            cd12 = self.wcs.cd[0][1]
            cd22 = self.wcs.cd[1][1]
            self.orientat = np.rad2deg(np.arctan2(cd12,cd22))
        except AttributeError:
            if self.wcs.has_cd():
                print "This file has a PC matrix. You may want to convert it \n \
                to a CD matrix, if reasonable, by running pc2.cd() method.\n \
                The orientation can be set then by calling setOrient() method.\n"
            self.pscale = None

    def updatePscale(self, scale):
        """
        Updates the CD matrix with a new plate scale
        """
        self.wcs.cd = self.wcs.cd/self.pscale*scale
        self.setPscale()

    def readModel(self, update=False, header=None):
        """
        Reads distortion model from IDCTAB.

        If IDCTAB is not found ('N/A', "", or not found on disk), then
        if SIP coefficients and first order IDCTAB coefficients are present
        in the header, restore the idcmodel from the header.
        If not - assign None to self.idcmodel.

        Parameters
        ----------
        header: pyfits.Header
                fits extension header
        update: boolean (False)
                if True - record the following IDCTAB quantities as header keywords:
                CX10, CX11, CY10, CY11, IDCSCALE, IDCTHETA, IDCXREF, IDCYREF,
                IDCV2REF, IDCV3REF
        """
        if self.idctab in [None, '', ' ','N/A']:
            #Keyword idctab is not present in header - check for sip coefficients
            if header is not None and 'IDCSCALE' in header:
                self._readModelFromHeader(header)
            else:
                print "Distortion model is not available: IDCTAB=None\n"
                self.idcmodel = None
        elif not os.path.exists(fileutil.osfn(self.idctab)):
            if header is not None and 'IDCSCALE' in header:
                self._readModelFromHeader(header)
            else:
                print 'Distortion model is not available: IDCTAB file %s not found\n' % self.idctab
                self.idcmodel = None
        else:
            self.readModelFromIDCTAB(header=header, update=update)

    def _readModelFromHeader(self, header):
        # Recreate idc model from SIP coefficients and header kw
        print 'Restoring IDC model from SIP coefficients\n'
        model = models.GeometryModel()
        cx, cy = coeff_converter.sip2idc(self)
        model.cx = cx
        model.cy = cy
        model.name = "sip"
        model.norder = header['A_ORDER']

        refpix = {}
        refpix['XREF'] = header['IDCXREF']
        refpix['YREF'] = header['IDCYREF']
        refpix['PSCALE'] = header['IDCSCALE']
        refpix['V2REF'] = header['IDCV2REF']
        refpix['V3REF'] = header['IDCV3REF']
        refpix['THETA'] = header['IDCTHETA']
        model.refpix = refpix

        self.idcmodel = model


    def readModelFromIDCTAB(self, header=None, update=False):
        """
        Read distortion model from idc table.

        Parameters
        ----------
        header: pyfits.Header
                fits extension header
        update: boolean (False)
                if True - save teh following as header keywords:
                CX10, CX11, CY10, CY11, IDCSCALE, IDCTHETA, IDCXREF, IDCYREF,
                IDCV2REF, IDCV3REF

        """
        if self.date_obs == None:
            print 'date_obs not available\n'
            self.idcmodel = None
            return
        if self.filter1 == None and self.filter2 == None:
            'No filter information available\n'
            self.idcmodel = None
            return

        self.idcmodel = models.IDCModel(self.idctab,
                    chip=self.chip, direction='forward', date=self.date_obs,
                    filter1=self.filter1, filter2=self.filter2,
                    offtab=self.offtab, binned=self.binned)

        if self.ltv1 != 0. or self.ltv2 != 0.:
            self.resetLTV()

        if update:
            if header==None:
                print 'Update header with IDC model kw requested but header was not provided\n.'
            else:
                self._updatehdr(header)

    def resetLTV(self):
        """
        Reset LTV values for polarizer data

        The polarizer field is smaller than the detector field.
        The distortion coefficients are defined for the entire
        polarizer field and the LTV values are set as with subarray
        data. This may also be true for other special filters.
        This is a special case when the observation is considered
        a subarray in terms of detector field but a full frame in
        terms of distortion model.
        To avoid shifting the distortion coefficients the LTV values
        are reset to 0.
        """
        if self.naxis1 == self.idcmodel.refpix['XSIZE'] and  \
           self.naxis2 == self.idcmodel.refpix['YSIZE']:
            self.ltv1 = 0.
            self.ltv2 = 0.

    def wcs2header(self, sip2hdr=False, idc2hdr=True, wcskey=None, relax=False):
        """
        Create a pyfits.Header object from WCS keywords.

        If the original header had a CD matrix, return a CD matrix,
        otherwise return a PC matrix.

        Parameters
        ----------
        sip2hdr: boolean
                 If True - include SIP coefficients
        """

        h = self.to_header(wkey=wcskey, relax=relax)
        if not wcskey:
            wcskey = self.wcs.alt
        if self.wcs.has_cd():
            h = altwcs.pc2cd(h, key=wcskey)

        if 'wcsname' not in h:
            if self.idctab is not None:
                wname = build_default_wcsname(self.idctab)
            else:
                wname = 'DEFAULT'
            h.update('wcsname'+wcskey, value=wname)

        if idc2hdr:
            for card in self._idc2hdr():
                h.update(card.key+wcskey, value=card.value, comment=card.comment)
        try:
            del h['RESTFRQ']
            del h['RESTWAV']
        except KeyError: pass

        if sip2hdr and self.sip:
            for card in self._sip2hdr('a'):
                h.update(card.key,value=card.value,comment=card.comment)
            for card in self._sip2hdr('b'):
                h.update(card.key,value=card.value,comment=card.comment)

            try:
                ap = self.sip.ap
            except AssertionError:
                ap = None
            try:
                bp = self.sip.bp
            except AssertionError:
                bp = None

            if ap:
                for card in self._sip2hdr('ap'):
                    h.update(card.key,value=card.value,comment=card.comment)
            if bp:
                for card in self._sip2hdr('bp'):
                    h.update(card.key,value=card.value,comment=card.comment)
        return h

    def _sip2hdr(self, k):
        """
        Get a set of SIP coefficients in the form of an array
        and turn them into a pyfits.Cardlist.
        k - one of 'a', 'b', 'ap', 'bp'
        """

        cards = pyfits.CardList()
        korder = self.sip.__getattribute__(k+'_order')
        cards.append(pyfits.Card(key=k.upper()+'_ORDER', value=korder))
        coeffs = self.sip.__getattribute__(k)
        ind = coeffs.nonzero()
        for i in range(len(ind[0])):
            card = pyfits.Card(key=k.upper()+'_'+str(ind[0][i])+'_'+str(ind[1][i]),
                               value=coeffs[ind[0][i], ind[1][i]])
            cards.append(card)
        return cards

    def _idc2hdr(self):
        # save some of the idc coefficients
        coeffs = ['ocx10', 'ocx11', 'ocy10', 'ocy11', 'idcscale']
        cards = pyfits.CardList()
        for c in coeffs:
            try:
                val = self.__getattribute__(c)
            except AttributeError:
                continue
            if val:
                cards.append(pyfits.Card(key=c, value=val))
        return cards

    def pc2cd(self):
        self.wcs.cd = self.wcs.pc.copy()

    def all_sky2pix(self,*args):
        """
        Performs full inverse transformation using iterative solution
        on full forward transformation with complete distortion model.

        NOTES
        -----
        Inputs can either be (RA, Dec, origin) or (RADec, origin) where RA and Dec
        are 1-D arrays/lists of coordinates and RADec is an array/list of pairs
        of coordinates.

        We now need to find the position we want by iterative
        improvement of an initial guess - the centre of the chip

        The method is to derive an "effective CD matrix" and use that
        to apply a correction until we are close enough (as defined by
        the ERR variable)

        Code from the drizzle task TRANBACK (dither$drizzle/tranback.f)
        defined the algorithm for this implementation

        """
        from stwcs.distortion import utils

        if len(args) == 2:
            xy, origin = args
            try:
                xy = np.asarray(xy)
                ra = xy[:,0]
                dec = xy[:,1]
                origin = int(origin)
            except:
                raise TypeError(
                    "When providing two arguments, they must be (RADec, origin)")
        elif len(args) == 3:
            ra, dec, origin = args
            try:
                ra = np.asarray(ra)
                dec = np.asarray(dec)
                origin = int(origin)
            except:
                raise TypeError(
                    "When providing three arguments, they must be (RA, Dec, origin)")
            if ra.size != dec.size:
                raise ValueError("RA and Dec arrays are not the same size")
        else:
            raise TypeError("Expected 2 or 3 arguments, %d given" % len(args))

        # Define some output arrays
        xout = np.zeros(len(ra),dtype=np.float64)
        yout = np.zeros(len(ra),dtype=np.float64)
        # ... and internal arrays
        x = np.zeros(3,dtype=np.float64)
        y = np.zeros(3,dtype=np.float64)

        # define delta for comparison
        err = 0.0001

        # Use linear WCS as frame in which to perform fit
        # rather than on the sky
        undistort = True
        if self.sip is None:
            # Only apply distortion if distortion coeffs are present.
            undistort = False
        wcslin = utils.output_wcs([self],undistort=undistort)

        # We can only transform 1 position at a time
        for r,d,n in zip(ra,dec,xrange(len(ra))):
            # First guess for output
            x[0],y[0] = self.wcs_sky2pix(r,d,origin)
            # also convert RA,Dec into undistorted linear pixel positions
            lx,ly = wcslin.wcs_sky2pix(r,d,origin)

            # Loop around until we are close enough (max 20 iterations)
            ev_old = None
            for i in xrange(20):
                x[1] = x[0] + 1.0
                y[1] = y[0]
                x[2] = x[0]
                y[2] = y[0] + 1.0
                # Perform full transformation on pixel position
                rao,deco = self.all_pix2sky(x,y,origin)
                # convert RA,Dec into linear pixel positions for fitting
                xo,yo = wcslin.wcs_sky2pix(rao,deco,origin)

                # Compute deltas between output and initial guess as
                # an affine transform then invert that transformation
                dxymat = np.array([[xo[1]-xo[0],yo[1]-yo[0]],
                          [xo[2]-xo[0],yo[2]-yo[0]]],dtype=np.float64)

                #invmat = np.linalg.inv(dxymat)
                # compute error in output position
                dx = lx - xo[0]
                dy = ly - yo[0]

                # record the old position
                xold = x[0]
                yold = y[0]

                # Update the initial guess position using the transform
                x[0] = xold + dx*dxymat[0][0] + dy*dxymat[1][0]
                y[0] = yold + dx*dxymat[0][1] + dy*dxymat[1][1]

                # Work out the error vector length
                ev = np.sqrt((x[0] - xold)**2 + (y[0] - yold)**2)

                # initialize record of previous error measurement during 1st iteration
                if ev_old is None:
                    ev_old = ev

                # Check to see whether we have reached the limit or
                # the new error is greater than error from previous iteration
                if ev < err or (np.abs(ev) > np.abs(ev_old)):
                    x[0] = xold
                    y[0] = yold
                    break
                # remember error measurement from previous iteration
                ev_old = ev

            xout[n] = x[0]
            yout[n] = y[0]

        return xout,yout

    def _updatehdr(self, ext_hdr):
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
        """
        Print the basic WCS keywords.
        """
        print 'WCS Keywords\n'
        print 'CD_11  CD_12: %r %r' % (self.wcs.cd[0,0],  self.wcs.cd[0,1])
        print 'CD_21  CD_22: %r %r' % (self.wcs.cd[1,0],  self.wcs.cd[1,1])
        print 'CRVAL    : %r %r' % (self.wcs.crval[0], self.wcs.crval[1])
        print 'CRPIX    : %r %r' % (self.wcs.crpix[0], self.wcs.crpix[1])
        print 'NAXIS    : %d %d' % (self.naxis1, self.naxis2)
        print 'Plate Scale : %r' % self.pscale
        print 'ORIENTAT : %r' % self.orientat
