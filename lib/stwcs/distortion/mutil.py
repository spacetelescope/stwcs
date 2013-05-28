from __future__ import division # confidence high

from stsci.tools import fileutil
import numpy as np
import string
import calendar

# Set up IRAF-compatible Boolean values
yes = True
no = False

# This function read the IDC table and generates the two matrices with
# the geometric correction coefficients.
#
#       INPUT: FITS object of open IDC table
#       OUTPUT: coefficient matrices for Fx and Fy
#
#### If 'tabname' == None: This should return a default, undistorted
####                        solution.
#

def readIDCtab (tabname, chip=1, date=None, direction='forward',
                filter1=None,filter2=None, offtab=None):

    """
        Read IDCTAB, and optional OFFTAB if sepcified, and generate
        the two matrices with the geometric correction coefficients.

        If tabname == None, then return a default, undistorted solution.
        If offtab is specified, dateobs also needs to be given.

    """

 # Return a default geometry model if no IDCTAB filename
    # is given.  This model will not distort the data in any way.
    if tabname == None:
        print 'Warning: No IDCTAB specified! No distortion correction will be applied.'
        return defaultModel()

    # Implement default values for filters here to avoid the default
    # being overwritten by values of None passed by user.
    if filter1 == None or filter1.find('CLEAR') == 0 or filter1.strip() == '':
        filter1 = 'CLEAR'
    if filter2 == None or filter2.find('CLEAR') == 0 or filter2.strip() == '':
        filter2 = 'CLEAR'

    # Insure that tabname is full filename with fully expanded
    # IRAF variables; i.e. 'jref$mc41442gj_idc.fits' should get
    # expanded to '/data/cdbs7/jref/mc41442gj_idc.fits' before
    # being used here.
    # Open up IDC table now...
    try:
        ftab = fileutil.openImage(tabname)
    except:
        err_str =  "------------------------------------------------------------------------ \n"
        err_str += "WARNING: the IDCTAB geometric distortion file specified in the image     \n"
        err_str += "header was not found on disk. Please verify that your environment        \n"
        err_str += "variable ('jref'/'uref'/'oref'/'nref') has been correctly defined. If    \n"
        err_str += "you do not have the IDCTAB file, you may obtain the latest version       \n"
        err_str += "of it from the relevant instrument page on the STScI HST website:        \n"
        err_str += "http://www.stsci.edu/hst/ For WFPC2, STIS and NICMOS data, the           \n"
        err_str += "present run will continue using the old coefficients provided in         \n"
        err_str += "the Dither Package (ca. 1995-1998).                                      \n"
        err_str += "------------------------------------------------------------------------ \n"
        raise IOError,err_str

    #First thing we need, is to read in the coefficients from the IDC
    # table and populate the Fx and Fy matrices.

    if 'DETECTOR' in ftab['PRIMARY'].header:
        detector = ftab['PRIMARY'].header['DETECTOR']
    else:
        if 'CAMERA' in ftab['PRIMARY'].header:
            detector = str(ftab['PRIMARY'].header['CAMERA'])
        else:
            detector = 1
    # First, read in TDD coeffs if present
    phdr = ftab['PRIMARY'].header
    instrument = phdr['INSTRUME']
    if instrument == 'ACS' and detector == 'WFC':
        skew_coeffs = read_tdd_coeffs(phdr)
    else:
        skew_coeffs = None

    # Set default filters for SBC
    if detector == 'SBC':
        if filter1 == 'CLEAR':
            filter1 = 'F115LP'
            filter2 = 'N/A'
        if filter2 == 'CLEAR':
            filter2 = 'N/A'

    # Read FITS header to determine order of fit, i.e. k
    norder = ftab['PRIMARY'].header['NORDER']
    if norder < 3:
        order = 3
    else:
        order = norder

    fx = np.zeros(shape=(order+1,order+1),dtype=np.float64)
    fy = np.zeros(shape=(order+1,order+1),dtype=np.float64)

    #Determine row from which to get the coefficients.
    # How many rows do we have in the table...
    fshape = ftab[1].data.shape
    colnames = ftab[1].data.names
    row = -1

    # Loop over all the rows looking for the one which corresponds
    # to the value of CCDCHIP we are working on...
    for i in xrange(fshape[0]):

        try:
            # Match FILTER combo to appropriate row,
            #if there is a filter column in the IDCTAB...
            if 'FILTER1' in colnames and 'FILTER2' in colnames:

                filt1 = ftab[1].data.field('FILTER1')[i]
                if filt1.find('CLEAR') > -1: filt1 = filt1[:5]

                filt2 = ftab[1].data.field('FILTER2')[i]
                if filt2.find('CLEAR') > -1: filt2 = filt2[:5]
            else:
                if 'OPT_ELEM' in colnames:
                    filt1 = ftab[1].data.field('OPT_ELEM')
                    if filt1.find('CLEAR') > -1: filt1 = filt1[:5]
                else:
                    filt1 = filter1

                if 'FILTER' in colnames:
                    _filt = ftab[1].data.field('FILTER')[i]
                    if _filt.find('CLEAR') > -1: _filt = _filt[:5]
                    if 'OPT_ELEM' in colnames:
                        filt2 = _filt
                    else:
                        filt1 = _filt
                        filt2 = 'CLEAR'
                else:
                    filt2 = filter2
        except:
            # Otherwise assume all rows apply and compare to input filters...
            filt1 = filter1
            filt2 = filter2

        if 'DETCHIP' in colnames:
            detchip = ftab[1].data.field('DETCHIP')[i]
            if not str(detchip).isdigit():
                detchip = 1
        else:
            detchip = 1

        if 'DIRECTION' in colnames:
            direct = string.strip(string.lower(ftab[1].data.field('DIRECTION')[i]))
        else:
            direct = 'forward'

        if filt1 == filter1.strip() and filt2 == filter2.strip():
            if direct == direction.strip():
                if int(detchip) == int(chip) or int(detchip) == -999:
                    row = i
                    break

    joinstr = ','
    if 'CLEAR' in filter1:
        f1str = ''
        joinstr = ''
    else:
        f1str = filter1.strip()
    if 'CLEAR' in filter2:
        f2str = ''
        joinstr = ''
    else:
        f2str = filter2.strip()
    filtstr = (joinstr.join([f1str,f2str])).strip()
    if row < 0:
        err_str = '\nProblem finding row in IDCTAB! Could not find row matching:\n'
        err_str += '        CHIP: '+str(detchip)+'\n'
        err_str += '     FILTERS: '+filtstr+'\n'
        ftab.close()
        del ftab
        raise LookupError,err_str
    else:
        print '- IDCTAB: Distortion model from row',str(row+1),'for chip',detchip,':',filtstr

    # Read in V2REF and V3REF: this can either come from current table,
    # or from an OFFTAB if time-dependent (i.e., for WFPC2)
    theta = None
    if 'V2REF' in colnames:
        v2ref = ftab[1].data.field('V2REF')[row]
        v3ref = ftab[1].data.field('V3REF')[row]
    else:
        # Read V2REF/V3REF from offset table (OFFTAB)
        if offtab:
            v2ref,v3ref,theta = readOfftab(offtab, date, chip=detchip)
        else:
            v2ref = 0.0
            v3ref = 0.0

    if theta == None:
        if 'THETA' in colnames:
            theta = ftab[1].data.field('THETA')[row]
        else:
            theta = 0.0

    refpix = {}
    refpix['XREF'] = ftab[1].data.field('XREF')[row]
    refpix['YREF'] = ftab[1].data.field('YREF')[row]
    refpix['XSIZE'] = ftab[1].data.field('XSIZE')[row]
    refpix['YSIZE'] = ftab[1].data.field('YSIZE')[row]
    refpix['PSCALE'] = round(ftab[1].data.field('SCALE')[row],8)
    refpix['V2REF'] = v2ref
    refpix['V3REF'] = v3ref
    refpix['THETA'] = theta
    refpix['XDELTA'] = 0.0
    refpix['YDELTA'] = 0.0
    refpix['DEFAULT_SCALE'] = yes
    refpix['centered'] = no
    refpix['skew_coeffs'] = skew_coeffs

    # Now that we know which row to look at, read coefficients into the
    #   numeric arrays we have set up...
    # Setup which column name convention the IDCTAB follows
    # either: A,B or CX,CY
    if 'CX10' in ftab[1].data.names:
        cxstr = 'CX'
        cystr = 'CY'
    else:
        cxstr = 'A'
        cystr = 'B'

    for i in xrange(norder+1):
        if i > 0:
            for j in xrange(i+1):
                xcname = cxstr+str(i)+str(j)
                ycname = cystr+str(i)+str(j)
                fx[i,j] = ftab[1].data.field(xcname)[row]
                fy[i,j] = ftab[1].data.field(ycname)[row]

    ftab.close()
    del ftab

    # If CX11 is 1.0 and not equal to the PSCALE, then the
    # coeffs need to be scaled

    if fx[1,1] == 1.0 and abs(fx[1,1]) != refpix['PSCALE']:
        fx *= refpix['PSCALE']
        fy *= refpix['PSCALE']

    # Return arrays and polynomial order read in from table.
    # NOTE: XREF and YREF are stored in Fx,Fy arrays respectively.
    return fx,fy,refpix,order
#
#
# Time-dependent skew correction coefficients (only ACS/WFC)
#
#
def read_tdd_coeffs(phdr):
    ''' Read in the TDD related keywords from the PRIMARY header of the IDCTAB
    '''
    skew_coeffs = {}
    skew_coeffs['TDDORDER'] = 0
    skew_coeffs['TDD_DATE'] = ""
    skew_coeffs['TDD_A'] = None
    skew_coeffs['TDD_B'] = None

    if "TDDORDER" in phdr:
        n = int(phdr["TDDORDER"])
    else:
        print 'TDDORDER kw not present, using default TDD correction'
        return None

    a = np.zeros((n+1,), np.float64)
    b = np.zeros((n+1,), np.float64)
    for i in range(n+1):
        a[i] = phdr.get(("TDD_A%d" % i), 0.0)
        b[i] = phdr.get(("TDD_B%d" % i), 0.0)
    if (a==0).all() and (b==0).all():
        print 'Warning: TDD_A and TDD_B coeffiecients have values of 0, \n \
               but TDDORDER is %d.' % TDDORDER

    skew_coeffs['TDDORDER'] = n
    skew_coeffs['TDD_DATE'] = phdr['TDD_DATE']
    skew_coeffs['TDD_A'] = a
    skew_coeffs['TDD_B'] = b

    return skew_coeffs

def readOfftab(offtab, date, chip=None):


#Read V2REF,V3REF from a specified offset table (OFFTAB).
# Return a default geometry model if no IDCTAB filenam  e
# is given.  This model will not distort the data in any way.

    if offtab == None:
        return 0.,0.

    # Provide a default value for chip
    if chip:
        detchip = chip
    else:
        detchip = 1

    # Open up IDC table now...
    try:
        ftab = fileutil.openImage(offtab)
    except:
        raise IOError,"Offset table '%s' not valid as specified!" % offtab

    #Determine row from which to get the coefficients.
    # How many rows do we have in the table...
    fshape = ftab[1].data.shape
    colnames = ftab[1].data.names
    row = -1

    row_start = None
    row_end = None

    v2end = None
    v3end = None
    date_end = None
    theta_end = None

    num_date = convertDate(date)
    # Loop over all the rows looking for the one which corresponds
    # to the value of CCDCHIP we are working on...
    for ri in xrange(fshape[0]):
        i = fshape[0] - ri - 1
        if 'DETCHIP' in colnames:
            detchip = ftab[1].data.field('DETCHIP')[i]
        else:
            detchip = 1

        obsdate = convertDate(ftab[1].data.field('OBSDATE')[i])

        # If the row is appropriate for the chip...
            # Interpolate between dates
        if int(detchip) == int(chip) or int(detchip) == -999:
            if num_date <= obsdate:
                date_end = obsdate
                v2end = ftab[1].data.field('V2REF')[i]
                v3end = ftab[1].data.field('V3REF')[i]
                theta_end = ftab[1].data.field('THETA')[i]
                row_end = i
                continue

            if row_end == None and (num_date > obsdate):
                date_end = obsdate
                v2end = ftab[1].data.field('V2REF')[i]
                v3end = ftab[1].data.field('V3REF')[i]
                theta_end = ftab[1].data.field('THETA')[i]
                row_end = i
                continue

            if num_date > obsdate:
                date_start = obsdate
                v2start = ftab[1].data.field('V2REF')[i]
                v3start = ftab[1].data.field('V3REF')[i]
                theta_start = ftab[1].data.field('THETA')[i]
                row_start = i
                break

    ftab.close()
    del ftab

    if row_start == None and row_end == None:
        print 'Row corresponding to DETCHIP of ',detchip,' was not found!'
        raise LookupError
    elif row_start == None:
        print '- OFFTAB: Offset defined by row',str(row_end+1)
    else:
        print '- OFFTAB: Offset interpolated from rows',str(row_start+1),'and',str(row_end+1)

    # Now, do the interpolation for v2ref, v3ref, and theta
    if row_start == None or row_end == row_start:
        # We are processing an observation taken after the last calibration
        date_start = date_end
        v2start = v2end
        v3start = v3end
        _fraction = 0.
        theta_start = theta_end
    else:
        _fraction = float((num_date - date_start)) / float((date_end - date_start))

    v2ref = _fraction * (v2end - v2start) + v2start
    v3ref = _fraction * (v3end - v3start) + v3start
    theta = _fraction * (theta_end - theta_start) + theta_start

    return v2ref,v3ref,theta

def readWCSCoeffs(header):

    #Read distortion coeffs from WCS header keywords and
    #populate distortion coeffs arrays.

    # Read in order for polynomials
    _xorder = header['a_order']
    _yorder = header['b_order']
    order = max(max(_xorder,_yorder),3)

    fx = np.zeros(shape=(order+1,order+1),dtype=np.float64)
    fy = np.zeros(shape=(order+1,order+1),dtype=np.float64)

    # Read in CD matrix
    _cd11 = header['cd1_1']
    _cd12 = header['cd1_2']
    _cd21 = header['cd2_1']
    _cd22 = header['cd2_2']
    _cdmat = np.array([[_cd11,_cd12],[_cd21,_cd22]])
    _theta = np.arctan2(-_cd12,_cd22)
    _rotmat = np.array([[np.cos(_theta),np.sin(_theta)],
                      [-np.sin(_theta),np.cos(_theta)]])
    _rCD = np.dot(_rotmat,_cdmat)
    _skew = np.arcsin(-_rCD[1][0] / _rCD[0][0])
    _scale = _rCD[0][0] * np.cos(_skew) * 3600.
    _scale2 = _rCD[1][1] * 3600.

    # Set up refpix
    refpix = {}
    refpix['XREF'] = header['crpix1']
    refpix['YREF'] = header['crpix2']
    refpix['XSIZE'] = header['naxis1']
    refpix['YSIZE'] = header['naxis2']
    refpix['PSCALE'] = _scale
    refpix['V2REF'] = 0.
    refpix['V3REF'] = 0.
    refpix['THETA'] = np.rad2deg(_theta)
    refpix['XDELTA'] = 0.0
    refpix['YDELTA'] = 0.0
    refpix['DEFAULT_SCALE'] = yes
    refpix['centered'] = yes


    # Set up template for coeffs keyword names
    cxstr = 'A_'
    cystr = 'B_'
    # Read coeffs into their own matrix
    for i in xrange(_xorder+1):
        for j in xrange(i+1):
            xcname = cxstr+str(j)+'_'+str(i-j)
            if xcname in header:
                fx[i,j] = header[xcname]

    # Extract Y coeffs separately as a different order may
    # have been used to fit it.
    for i in xrange(_yorder+1):
        for j in xrange(i+1):
            ycname = cystr+str(j)+'_'+str(i-j)
            if ycname in header:
                fy[i,j] = header[ycname]

    # Now set the linear terms
    fx[0][0] = 1.0
    fy[0][0] = 1.0

    return fx,fy,refpix,order


def readTraugerTable(idcfile,wavelength):

    # Return a default geometry model if no coefficients filename
    # is given.  This model will not distort the data in any way.
    if idcfile == None:
        return fileutil.defaultModel()

    # Trauger coefficients only result in a cubic file...
    order = 3
    numco = 10
    a_coeffs = [0] * numco
    b_coeffs = [0] * numco
    indx = _MgF2(wavelength)

    ifile = open(idcfile,'r')
    # Search for the first line of the coefficients
    _line = fileutil.rAsciiLine(ifile)
    while string.lower(_line[:7]) != 'trauger':
        _line = fileutil.rAsciiLine(ifile)
    # Read in each row of coefficients,split them into their values,
    # and convert them into cubic coefficients based on
    # index of refraction value for the given wavelength
    # Build X coefficients from first 10 rows of Trauger coefficients
    j = 0
    while j < 20:
        _line = fileutil.rAsciiLine(ifile)
        if _line == '': continue
        _lc = string.split(_line)
        if j < 10:
            a_coeffs[j] = float(_lc[0])+float(_lc[1])*(indx-1.5)+float(_lc[2])*(indx-1.5)**2
        else:
            b_coeffs[j-10] = float(_lc[0])+float(_lc[1])*(indx-1.5)+float(_lc[2])*(indx-1.5)**2
        j = j + 1

    ifile.close()
    del ifile

    # Now, convert the coefficients into a Numeric array
    # with the right coefficients in the right place.
    # Populate output values now...
    fx = np.zeros(shape=(order+1,order+1),dtype=np.float64)
    fy = np.zeros(shape=(order+1,order+1),dtype=np.float64)
    # Assign the coefficients to their array positions
    fx[0,0] = 0.
    fx[1] = np.array([a_coeffs[2],a_coeffs[1],0.,0.],dtype=np.float64)
    fx[2] = np.array([a_coeffs[5],a_coeffs[4],a_coeffs[3],0.],dtype=np.float64)
    fx[3] = np.array([a_coeffs[9],a_coeffs[8],a_coeffs[7],a_coeffs[6]],dtype=np.float64)
    fy[0,0] = 0.
    fy[1] = np.array([b_coeffs[2],b_coeffs[1],0.,0.],dtype=np.float64)
    fy[2] = np.array([b_coeffs[5],b_coeffs[4],b_coeffs[3],0.],dtype=np.float64)
    fy[3] = np.array([b_coeffs[9],b_coeffs[8],b_coeffs[7],b_coeffs[6]],dtype=np.float64)

    # Used in Pattern.computeOffsets()
    refpix = {}
    refpix['XREF'] = None
    refpix['YREF'] = None
    refpix['V2REF'] = None
    refpix['V3REF'] = None
    refpix['XDELTA'] = 0.
    refpix['YDELTA'] = 0.
    refpix['PSCALE'] = None
    refpix['DEFAULT_SCALE'] = no
    refpix['centered'] = yes

    return fx,fy,refpix,order


def readCubicTable(idcfile):
    # Assumption: this will only be used for cubic file...
    order = 3
    # Also, this function does NOT perform any scaling on
    # the coefficients, it simply passes along what is found
    # in the file as is...

    # Return a default geometry model if no coefficients filename
    # is given.  This model will not distort the data in any way.
    if idcfile == None:
        return fileutil.defaultModel()

    ifile = open(idcfile,'r')
    # Search for the first line of the coefficients
    _line = fileutil.rAsciiLine(ifile)

    _found = no
    while _found == no:
        if _line[:7] in  ['cubic','quartic','quintic'] or _line[:4] == 'poly':
            found = yes
            break
        _line = fileutil.rAsciiLine(ifile)

    # Read in each row of coefficients, without line breaks or newlines
    # split them into their values, and create a list for A coefficients
    # and another list for the B coefficients
    _line = fileutil.rAsciiLine(ifile)
    a_coeffs = string.split(_line)

    x0 = float(a_coeffs[0])
    _line = fileutil.rAsciiLine(ifile)
    a_coeffs[len(a_coeffs):] = string.split(_line)
    # Scale coefficients for use within PyDrizzle
    for i in range(len(a_coeffs)):
        a_coeffs[i] = float(a_coeffs[i])

    _line = fileutil.rAsciiLine(ifile)
    b_coeffs = string.split(_line)
    y0 = float(b_coeffs[0])
    _line = fileutil.rAsciiLine(ifile)
    b_coeffs[len(b_coeffs):] = string.split(_line)
    # Scale coefficients for use within PyDrizzle
    for i in range(len(b_coeffs)):
        b_coeffs[i] = float(b_coeffs[i])

    ifile.close()
    del ifile
    # Now, convert the coefficients into a Numeric array
    # with the right coefficients in the right place.
    # Populate output values now...
    fx = np.zeros(shape=(order+1,order+1),dtype=np.float64)
    fy = np.zeros(shape=(order+1,order+1),dtype=np.float64)
    # Assign the coefficients to their array positions
    fx[0,0] = 0.
    fx[1] = np.array([a_coeffs[2],a_coeffs[1],0.,0.],dtype=np.float64)
    fx[2] = np.array([a_coeffs[5],a_coeffs[4],a_coeffs[3],0.],dtype=np.float64)
    fx[3] = np.array([a_coeffs[9],a_coeffs[8],a_coeffs[7],a_coeffs[6]],dtype=np.float64)
    fy[0,0] = 0.
    fy[1] = np.array([b_coeffs[2],b_coeffs[1],0.,0.],dtype=np.float64)
    fy[2] = np.array([b_coeffs[5],b_coeffs[4],b_coeffs[3],0.],dtype=np.float64)
    fy[3] = np.array([b_coeffs[9],b_coeffs[8],b_coeffs[7],b_coeffs[6]],dtype=np.float64)

    # Used in Pattern.computeOffsets()
    refpix = {}
    refpix['XREF'] = None
    refpix['YREF'] = None
    refpix['V2REF'] = x0
    refpix['V3REF'] = y0
    refpix['XDELTA'] = 0.
    refpix['YDELTA'] = 0.
    refpix['PSCALE'] = None
    refpix['DEFAULT_SCALE'] = no
    refpix['centered'] = yes

    return fx,fy,refpix,order

def factorial(n):
    """ Compute a factorial for integer n. """
    m = 1
    for i in range(int(n)):
        m = m * (i+1)
    return m

def combin(j,n):
    """ Return the combinatorial factor for j in n."""
    return (factorial(j) / (factorial(n) * factorial( (j-n) ) ) )


def defaultModel():
    """ This function returns a default, non-distorting model
        that can be used with the data.
    """
    order = 3

    fx = np.zeros(shape=(order+1,order+1),dtype=np.float64)
    fy = np.zeros(shape=(order+1,order+1),dtype=np.float64)

    fx[1,1] = 1.
    fy[1,0] = 1.

    # Used in Pattern.computeOffsets()
    refpix = {}
    refpix['empty_model'] = yes
    refpix['XREF'] = None
    refpix['YREF'] = None
    refpix['V2REF'] = 0.
    refpix['XSIZE'] = 0.
    refpix['YSIZE'] = 0.
    refpix['V3REF'] = 0.
    refpix['XDELTA'] = 0.
    refpix['YDELTA'] = 0.
    refpix['PSCALE'] = None
    refpix['DEFAULT_SCALE'] = no
    refpix['THETA'] = 0.
    refpix['centered'] = yes
    return fx,fy,refpix,order

# Function to compute the index of refraction for MgF2 at
# the specified wavelength for use with Trauger coefficients
def _MgF2(lam):
    _sig = pow((1.0e7/lam),2)
    return np.sqrt(1.0 + 2.590355e10/(5.312993e10-_sig) +
        4.4543708e9/(11.17083e9-_sig) + 4.0838897e5/(1.766361e5-_sig))


def convertDate(date):
    """ Converts the DATE-OBS date string into an integer of the
        number of seconds since 1970.0 using calendar.timegm().

        INPUT: DATE-OBS in format of 'YYYY-MM-DD'.
        OUTPUT: Date (integer) in seconds.
    """

    _dates = date.split('-')
    _val = 0
    _date_tuple = (int(_dates[0]), int(_dates[1]), int(_dates[2]), 0, 0, 0, 0, 0, 0)

    return calendar.timegm(_date_tuple)
