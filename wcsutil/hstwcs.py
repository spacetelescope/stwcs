from __future__ import division # confidence high

import os.path
from pywcs import WCS
import pyfits
import instruments
from stwcs.distortion import models, coeff_converter
import altwcs
import numpy as np
from pytools import fileutil
from pytools.fileutil import DEGTORAD, RADTODEG

import mappings
from mappings import inst_mappings, ins_spec_kw
from mappings import basic_wcs


__docformat__ = 'restructuredtext'
__version__ = '0.4'

class HSTWCS(WCS):
    """
    Purpose
    =======
    Create a WCS object based on the instrument.
    It has all basic WCS kw as attributes (set by pywcs).
    It also uses the primary and extension header to define
    instrument specific attributes.
    """
    def __init__(self, fobj='DEFAULT', ext=None, minerr=0.0, wcskey=" "):
        """
        :Parameters:
        `fobj`: string or PyFITS HDUList object or None
                a file name, e.g j9irw4b1q_flt.fits
                a fully qualified filename[EXTNAME,EXTNUM], e.g. j9irw4b1q_flt.fits[sci,1]
                a pyfits file object, e.g pyfits.open('j9irw4b1q_flt.fits'), in which case the 
                     user is responsible for closing the file object.
        `ext`:  int or None
                extension number
                if ext==None, it is assumed the data is in the primary hdu
 
        `minerr`: float
                minimum value a distortion correction must have in order to be applied.
                If CPERRja, CQERRja are smaller than minerr, the corersponding
                distortion is not applied.
        """

        self.inst_kw = ins_spec_kw
        self.minerr = minerr
        self.wcskey = wcskey
        
        if fobj != 'DEFAULT':
            filename, hdr0, ehdr, phdu = self.parseInput(f=fobj, ext=ext)
            self.filename = filename
            self.instrument = hdr0['INSTRUME']
            
            WCS.__init__(self, ehdr, fobj=phdu, minerr=self.minerr, key=self.wcskey)
            # If input was a pyfits HDUList object, it's the user's
            # responsibility to close it, otherwise, it's closed here.
            if not isinstance(fobj, pyfits.HDUList):
                phdu.close()   
            self.setInstrSpecKw(hdr0, ehdr)
            self.readIDCCoeffs(ehdr)
            extname = ehdr.get('EXTNAME', "")
            extnum = ehdr.get('EXTVER', None)
            self.extname = (extname, extnum)
        else:
            # create a default HSTWCS object 
            self.instrument = 'DEFAULT'
            WCS.__init__(self, minerr=self.minerr, key=self.wcskey)
            self.wcs.cd = np.array([[1.0, 0.0], [0.0, 1.0]], np.double)
            self.wcs.crval = np.zeros((self.naxis,), np.double)
            self.wcs.crpix = np.zeros((self.naxis,), np.double)
            self.wcs.set()
            self.setInstrSpecKw()
        self.setPscale()
        self.setOrient()
        
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
        self.wcs.set()


    def readModel(self, update=False, header=None):
        """
        Reads distortion model from IDCTAB.
        If IDCTAB is not found ('N/A', "", or not found on disk), then
        if SIP coefficients and first order IDCTAb coefficients are present
        in the header, restore the idcmodel from the header.
        If not - assign None to self.idcmodel.
        """

        if self.idctab == None or self.idctab == ' ':
            #Keyword idctab is not present in header - check for sip coefficients
            if header.has_key('IDCSCALE'):
                self._readModelFromHeader(header)
            else:
                print "Distortion model is not available: IDCTAB=None\n"
                self.idcmodel = None
        elif not os.path.exists(fileutil.osfn(self.idctab)):
            if header.has_key('IDCSCALE'):
                self._readModelFromHeader(header)
            else:
                print 'Distortion model is not available: IDCTAB file %s not found\n' % self.idctab
                self.idcmodel = None
        else:
            self._readModelFromIDCTAB(header=header, update=update)
            
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
        
        
    def _readModelFromIDCTAB(self, header=None, update=False):
        """
        Purpose
        =======
        Read distortion model from idc table.
        Save some of the information as kw needed for interpreting the distortion
        If header is provided and update is True, some IDC model kw
        will be recorded in the header.
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

        if update:
            if header==None:
                print 'Update header with IDC model kw requested but header was not provided\n.'
            else:
                self._updatehdr(header)

    def wcs2header(self, sip2hdr=False):
        """
        Create a pyfits.Header object from all WCS keywords,
        including the SIP coefficients.
        """
        h = self.to_header()
        if self.wcs.has_cd():
            h = altwcs.pc2cd(h)
            
        if sip2hdr:
            hwcs = h.ascardlist()
            cards = self._sip2hdr('a')
            hwcs.extend(cards)
            cards = self._sip2hdr('b')
            hwcs.extend(cards)
            
            try:
                ap = self.sip.ap
            except AssertionError:
                ap = None
            try:
                bp = self.sip.bp
            except AssertionError:
                bp = None
            
            if ap:
                cards = self._sip2hdr('ap')
                hwcs.extend(cards)
            if bp:
                cards = self._sip2hdr('bp')
                hwcs.extend(cards)
                
            h = pyfits.Header(hwcs)
        
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
        print 'WCS Keywords\n'
        print 'CD_11  CD_12: %r %r' % (self.wcs.cd[0,0],  self.wcs.cd[0,1])
        print 'CD_21  CD_22: %r %r' % (self.wcs.cd[1,0],  self.wcs.cd[1,1])
        print 'CRVAL    : %r %r' % (self.wcs.crval[0], self.wcs.crval[1])
        print 'CRPIX    : %r %r' % (self.wcs.crpix[0], self.wcs.crpix[1])
        print 'NAXIS    : %d %d' % (self.naxis1, self.naxis2)
        print 'Plate Scale : %r' % self.pscale
        print 'ORIENTAT : %r' % self.orientat


