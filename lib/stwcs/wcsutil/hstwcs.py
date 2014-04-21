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


class NoConvergence(Exception):
    """
    An error class used to report non-convergence and/or divergence of
    numerical methods. It is used to report errors in the iterative solution
    used by the :py:meth:`~stwcs.hstwcs.HSTWCS.all_sky2pix`\ .

    Attributes
    ----------

    best_solution : numpy.array
        Best solution achieved by the method.

    accuracy : float
        Accuracy of the :py:attr:`best_solution`\ .

    niter : int
        Number of iterations performed by the numerical method to compute
        :py:attr:`best_solution`\ .

    divergent : None, numpy.array
        Indices of the points in :py:attr:`best_solution` array for which the
        solution appears to be divergent. If the solution does not diverge,
        `divergent` will be set to `None`.

    nonconvergent : None, numpy.array
        Indices of the points in :py:attr:`best_solution` array for which the
        solution failed to converge within the specified maximum number
        of iterations. If there are no non-converging poits (i.e., if
        the required accuracy has been achieved for all points) then
        `nonconvergent` will be set to `None`.

    """
    def __init__(self, *args, **kwargs):
        super(NoConvergence, self).__init__(*args)

        self.best_solution  = kwargs.pop('best_solution', None)
        self.accuracy       = kwargs.pop('accuracy', None)
        self.niter          = kwargs.pop('niter', None)
        self.divergent      = kwargs.pop('divergent', None)
        self.nonconvergent  = kwargs.pop('nonconvergent', None)


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

    def all_sky2pix(self, *args, **kwargs):
        """
        all_sky2pix(*arg, accuracy=1.0e-4, maxiter=20, adaptive=False, quiet=False)

        Performs full inverse transformation using iterative solution
        on full forward transformation with complete distortion model.

        Parameters
        ----------
        accuracy : float, optional (Default = 1.0e-4)
            Required accuracy of the solution.

        maxiter : int, optional (Default = 20)
            Maximum number of iterations allowed to reach the solution.

        adaptive : bool, optional (Default = False)
            Specifies whether to adaptively select only points that did not
            converge to a solution whithin the required accuracy for the
            next iteration. Default is recommended for HST instruments.

            .. note::
               The :py:meth:`all_sky2pix` uses a vectorized implementation of
               the method of consecutive approximations (see `Notes` section
               below) in which it iterates over *all* input poits *regardless*
               untill the required accuracy has been reached for *all* input
               points. In some cases it may be possible that *almost all* points
               have reached the required accuracy but there are only a few
               of input data points for which additional iterations may be
               needed. In this situation it may be advantageous to set
               `adaptive` = `True`\ in which case :py:meth:`all_sky2pix` will
               continue iterating *only* over the points that have not yet
               converged to the required accuracy. However, for the HST's
               ACS/WFC detector, which has the strongest distortions of all
               HST instruments, testing has shown that enabling this option
               would lead to a 10-30\% penalty in computational time.
               Therefore, for HST instruments, it is recommended to set
               `adaptive` = `False`\ .

        quiet : bool, optional (Default = False)
            Do not throw :py:class:`NoConvergence` exceptions when the method
            does not converge to a solution with the required accuracy
            within a specified number of maximum iterations set by `maxiter`
            parameter. Instead, simply return the found solution.

        Raises
        ------
        NoConvergence
            The method does not converge to a
            solution with the required accuracy within a specified number
            of maximum iterations set by the `maxiter` parameter.

        Notes
        -----
        Inputs can either be (RA, Dec, origin) or (RADec, origin) where RA and Dec
        are 1-D arrays/lists of coordinates and RADec is an array/list of pairs
        of coordinates.

        Using the method of consecutive approximations we iterate starting
        with the initial approximation, which is computed using the
        non-distorion-aware :py:meth:`wcs_sky2pix` (or equivalent).

        The :py:meth:`all_sky2pix` function uses a vectorized implementation
        of the method of consecutive approximations and therefore it is
        highly efficient (>30x) when *all* data points that need to be
        converted from sky coordinates to image coordinates are passed at
        *once*\ . Therefore, it is advisable, whenever possible, to pass
        as input a long array of all points that need to be converted to
        :py:meth:`all_sky2pix` instead of calling :py:meth:`all_sky2pix`
        for each data point. Also see the note to the `adaptive` parameter.

        Examples
        --------
        >>> import stwcs, pyfits
        >>> hdulist = pyfits.open('j94f05bgq_flt.fits')
        >>> w = stwcs.wcsutil.HSTWCS(hdulist, ext=('sci',1))
        >>> hdulist.close()
        >>> ra, dec = w.all_pix2sky([1,2,3],[1,1,1],1); print(ra); print(dec)
        [ 5.52645241  5.52649277  5.52653313]
        [-72.05171776 -72.05171295 -72.05170814]
        >>> radec = w.all_pix2sky([[1,1],[2,1],[3,1]],1); print(radec)
        [[  5.52645241 -72.05171776]
         [  5.52649277 -72.05171295]
         [  5.52653313 -72.05170814]]
        >>> x, y = w.all_sky2pix(ra,dec,1)
        >>> print(x)
        [ 1.00000233  2.00000232  3.00000233]
        >>> print(y)
        [ 0.99999997  0.99999997  0.99999998]
        >>> xy = w.all_sky2pix(radec,1)
        >>> print(xy)
        [[ 1.00000233  0.99999997]
         [ 2.00000232  0.99999997]
         [ 3.00000233  0.99999998]]
        >>> xy = w.all_sky2pix(radec,1, maxiter=3, accuracy=1.0e-10, quiet=False)
        NoConvergence: 'HSTWCS.all_sky2pix' failed to converge to requested \
accuracy after 3 iterations.

        """
        nargs = len(args)

        if nargs == 3:
            try:
                ra     = np.asarray(args[0], dtype=np.float64)
                dec    = np.asarray(args[1], dtype=np.float64)
                #assert( len(ra.shape) == 1 and len(dec.shape) == 1 )
                origin = int(args[2])
                vect1D = True
            except:
                raise TypeError("When providing three arguments, they must " \
                            "be (Ra, Dec, origin) where Ra and Dec are " \
                            "Nx1 vectors.")
        elif nargs == 2:
            try:
                rd  = np.asarray(args[0], dtype=np.float64)
                #assert( rd.shape[1] == 2 )
                ra  = rd[:,0]
                dec = rd[:,1]
                origin = int(args[1])
                vect1D = False
            except:
                raise TypeError("When providing two arguments, they must " \
                            "be (RaDec, origin) where RaDec is a Nx2 array.")
        else:
            raise TypeError("Expected 2 or 3 arguments, {:d} given." \
                            .format(nargs))

        # process optional arguments:
        accuracy = kwargs.pop('accuracy', 1.0e-4)
        maxiter  = kwargs.pop('maxiter', 20)
        quiet    = kwargs.pop('quiet', False)
        adaptive = kwargs.pop('adaptive', False)

        # initialize iterative process:
        x0, y0 = self.wcs_sky2pix(ra, dec, origin) # initial approximation (WCS based only)

        # see if iterative solution is required (when any of the
        # non-CD-matrix corrections are present). If not required
        # return initial approximation (xy0).
        if self.sip is None and \
           self.cpdis1 is None and self.cpdis2 is None and \
           self.det2im1 is None and self.det2im2 is None:
            # no non-WCS corrections are detected - return initial approximation
            if vect1D:
                return [x0, y0]
            else:
                return np.dstack([x0,y0])[0]

        x  = x0.copy() # 0-order solution
        y  = y0.copy() # 0-order solution

        # initial correction:
        dx, dy = self.pix2foc(x, y, origin)
        # If pix2foc does not apply all the required distortion
        # corrections then replace the above line with:
        #r0, d0 = self.all_pix2sky(x, y, origin)
        #dx, dy = self.wcs_sky2pix(r0, d0, origin )
        dx -= x0
        dy -= y0

        # update initial solution:
        x -= dx
        y -= dy

        # norn (L2) squared of the correction:
        dn2prev = dx**2+dy**2
        dn2     = dn2prev

        # process all coordinates simultaneously:
        iterlist  = range(1, maxiter+1)
        accuracy2 = accuracy**2
        ind       = None
        inddiv    = None

        divergent = False

        if not adaptive:
            for k in iterlist:
                # check convergence:
                if np.max(dn2) < accuracy2:
                    break

                # check for divergence:
                inddiv, = np.where((dn2 > dn2prev) & (dn2 >= accuracy2))
                if inddiv.shape[0] > 0:
                    divergent = True
                    break

                # find correction to the previous solution:
                dx, dy = self.pix2foc(x, y, origin)
                # If pix2foc does not apply all the required distortion
                # corrections then replace the above line with:
                #r0, d0 = self.all_pix2sky(x, y, origin)
                #dx, dy = self.wcs_sky2pix(r0, d0, origin )

                dx -= x0
                dy -= y0

                # apply correction:
                x -= dx
                y -= dy

                # update norn (L2) squared of the correction:
                dn2prev = dn2.copy()
                dn2     = dx**2+dy**2

        else:
            ind, = np.where(dn2 >= accuracy2)

            for k in iterlist:
                # check convergence:
                if ind.shape[0] == 0:
                    break

                # check for divergence:
                inddiv = ind[np.where(dn2[ind] > dn2prev[ind])]
                if inddiv.shape[0] > 0:
                    divergent = True
                    break

                # find correction to the previous solution:
                dx[ind], dy[ind] = self.pix2foc(x[ind], y[ind], origin)
                # If pix2foc does not apply all the required distortion
                # corrections then replace the above line with:
                #r0[ind], d0[ind] = self.all_pix2sky(x[ind], y[ind], origin)
                #dx[ind], dy[ind] = self.wcs_sky2pix(r0[ind], d0[ind], origin )
                dx[ind] -= x0[ind]
                dy[ind] -= y0[ind]

                # apply correction:
                x[ind] -= dx[ind]
                y[ind] -= dy[ind]

                # update norn (L2) squared of the correction:
                dn2prev  = dn2.copy()
                dn2 = dx**2+dy**2

                # update indices of elements that still need correction:
                ind, = np.where(dn2 >= accuracy2)
                #ind = ind[np.where(dx[ind]**2+dy[ind]**2 >= accuracy2)]

        if k >= maxiter and not quiet:
            if vect1D:
                sol  = [x, y]
                err  = [np.abs(dx), np.abs(dy)]
            else:
                sol  = np.dstack( [x, y] )[0]
                err  = np.dstack( [np.abs(dx), np.abs(dy)] )[0]

            if ind is None:
                ind, = np.where(dx**2+dy**2 >= accuracy2)

            if inddiv is None:
                inddiv, = np.where(dn2[ind] > dn2prev[ind])

            if ind.shape[0] == 0:
                ind    = None
                inddiv = None

            elif inddiv.shape[0] == 0:
                inddiv = None

            assert(ind is not None or inddiv is not None) # <-- sanity check

            raise NoConvergence("'HSTWCS.all_sky2pix' failed to converge to "\
                                "requested accuracy after {:d} iterations."  \
                                .format(k), best_solution = sol,             \
                                accuracy = err, niter = k,                   \
                                nonconvergent = ind, divergent = inddiv)

        if vect1D:
            return [x, y]
        else:
            return np.dstack([x,y])[0]


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
