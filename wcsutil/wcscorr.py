import os,copy
import pyfits
import numpy as np

from pytools import fileutil
import stwcs
from stwcs.wcsutil import altwcs
import convertwcs


DEFAULT_WCS_KEYS = ['CRVAL1','CRVAL2','CRPIX1','CRPIX2',
                    'CD1_1','CD1_2','CD2_1','CD2_2',
                    'CTYPE1','CTYPE2']
DEFAULT_PRI_KEYS = ['PA_V3']
###
### WCSEXT table related keyword archive functions
###
def init_wcscorr(input, force=False):
    """
    This function will initialize the WCSCORR table if it is not already present,
    and look for WCS keywords with a prefix of 'O' as the original OPUS
    generated WCS as the initial row for the table or use the current WCS
    keywords as initial row if no 'O' prefix keywords are found.

    This function will NOT overwrite any rows already present.

    This function works on all SCI extensions at one time.
    """

    # TODO: Create some sort of decorator or (for Python2.5) context for
    # opening a FITS file and closing it when done, if necessary
    if not isinstance(input, pyfits.HDUList):
        # input must be a filename, so open as PyFITS object
        fimg = pyfits.open(input, mode='update')
        need_to_close = True
    else:
        fimg = input
        need_to_close = False

    # Verify that a WCSCORR extension does not already exist...
    for extn in fimg:
        if extn.header.has_key('extname') and \
           extn.header['extname'] == 'WCSCORR':
            if not force:
                return
            else:
                del fimg['WCSCORR']
    # define the primary columns of the WCSEXT table with initial rows for each
    # SCI extension for the original OPUS solution
    numsci = fileutil.countExtn(fimg)

    # create new table with more rows than needed initially to make it easier to
    # add new rows later
    wcsext = create_wcscorr(numrows=numsci, padding=numsci * 4)
    # Assign the correct EXTNAME value to this table extension
    wcsext.header.update('TROWS', numsci * 2,
                         comment='Number of updated rows in table')
    wcsext.header.update('EXTNAME', 'WCSCORR',
                         comment='Table with WCS Update history')

    used_wcskeys = None
    wcs1 = stwcs.wcsutil.HSTWCS(fimg,ext=('SCI',1))
    idc2header = True
    if wcs1.idcscale is None:
        idc2header = False
    wcs_keywords = wcs1.wcs2header(idc2hdr=idc2header).keys()

    # Now copy original OPUS values into table
    for extver in xrange(1, numsci + 1):
        rowind = find_wcscorr_row(wcsext.data,
                                  {'WCS_ID': 'OPUS', 'EXTVER': extver,
                                   'WCS_key':'O'})
        # There should only EVER be a single row for each extension with OPUS values
        rownum = np.where(rowind)[0][0]
        #print 'Archiving OPUS WCS in row number ',rownum,' in WCSCORR table for SCI,',extver

        hdr = fimg['SCI', extver].header
        # define set of WCS keywords which need to be managed and copied to the table
        if used_wcskeys is None:
            used_wcskeys = altwcs.wcskeys(hdr)
        # Check to see whether or not there is an OPUS alternate WCS present,
        # if so, get its values directly, otherwise, archive the PRIMARY WCS
        # as the OPUS values in the WCSCORR table
        if 'O' in used_wcskeys:
            wkey = 'O'
        else:
            wkey = ' '
        wcs = stwcs.wcsutil.HSTWCS(fimg, ext=('SCI', extver), wcskey=wkey)
        wcshdr = wcs.wcs2header(idc2hdr=idc2header)

        for key in DEFAULT_PRI_KEYS:
            prihdr_keys = []
            if not hdr.has_key(key):
                prihdr_keys.append(key)

        if wcsext.data.field('CRVAL1')[rownum] != 0:
            # If we find values for these keywords already in the table, do not
            # overwrite them again
            print 'WCS keywords already updated...'
            break
        for key in wcs_keywords:
            if key in wcsext.data.names:
                wcsext.data.field(key)[rownum] = wcshdr[(key+wkey)[:8]]
        # Now get any keywords from PRIMARY header needed for WCS updates
        for key in prihdr_keys:
            wcsext.data.field(key)[rownum] = fimg[0].header[key]

    # Now that we have archived the OPUS alternate WCS, remove it from the list
    # of used_wcskeys
    if 'O' in used_wcskeys:
        used_wcskeys.remove('O')

    # Now copy remaining alternate WCSs into table
    # TODO: Much of this appears to be redundant with update_wcscorr; consider
    # merging them...
    for uwkey in used_wcskeys:
        if wkey == ' ':
            break
        for extver in xrange(1, numsci + 1):
            hdr = fimg['SCI', extver].header
            wcs = stwcs.wcsutil.HSTWCS(fimg, ext=('SCI', extver),
                                       wcskey=uwkey)
            wcshdr = wcs.wcs2header()
            if 'WCSNAME' + uwkey not in wcshdr:
                idctab = fileutil.osfn(fimg[0].header['idctab'])
                idcname = os.path.split(idctab)[-1]
                idcname = idcname[:idcname.find('_')]
                wcsid = 'IDC_' + idcname
            else:
                wcsid = wcshdr['WCSNAME' + uwkey]

            # identify next empty row
            rowind = find_wcscorr_row(wcsext.data, selections={'wcs_id': ''})
            rows = np.where(rowind)
            if len(rows[0]) > 0:
                rownum = np.where(rowind)[0][0]
            else:
                print 'No available rows found for updating. '
            #print 'Archiving current WCS row number ',rownum,' in WCSCORR table for SCI,',extver

            # Update selection columns for this row with relevant values
            wcsext.data.field('WCS_ID')[rownum] = wcsid
            wcsext.data.field('EXTVER')[rownum] = extver
            wcsext.data.field('WCS_key')[rownum] = uwkey

            # Look for standard WCS keyword values
            for key in wcs_keywords:
                if key in wcsext.data.names:
                    wcsext.data.field(key)[rownum] = wcshdr[key + uwkey]
            # Now get any keywords from PRIMARY header needed for WCS updates
            for key in prihdr_keys:
                wcsext.data.field(key)[rownum] = fimg[0].header[key]

    # Append this table to the image FITS file
    fimg.append(wcsext)
    # force an update now
    # set the verify flag to 'warn' so that it will always succeed, but still
    # tell the user if PyFITS detects any problems with the file as a whole
    fimg.flush('warn')

    if need_to_close:
        fimg.close()


def find_wcscorr_row(wcstab, selections):
    """
    Return an array of indices from the table (NOT HDU) 'wcstab' that matches the
    selections specified by the user.

    The row selection criteria must be specified as a dictionary with
    column name as key and value(s) representing the valid desired row values.
    For example, {'wcs_id':'OPUS','extver':2}.
    """

    mask = None
    for i in selections:
        icol = wcstab.field(i)
        if isinstance(icol,np.chararray): icol = icol.rstrip()
        bmask = (icol == selections[i])
        if mask is None:
            mask = bmask.copy()
        else:
            mask = np.logical_and(mask,bmask)
        del bmask
    return mask


def archive_wcs_file(image, wcs_id=None):
    """
    Update WCSCORR table with rows for each SCI extension to record the
    newly updated WCS keyword values.
    """

    if not isinstance(image, pyfits.HDUList):
        fimg = pyfits.open(image, mode='update')
        close_image = True
    else:
        fimg = image
        close_image = False

    update_wcscorr(fimg, wcs_id=wcs_id)

    if close_image:
        fimg.close()


def update_wcscorr(dest, source=None, extname='SCI', wcs_id=None):
    """
    Update WCSCORR table with a new row or rows for this extension header. It
    copies the current set of WCS keywords as a new row of the table based on
    keyed WCSs as per Paper I Multiple WCS standard).

    Parameters
    ----------
    dest : HDUList
        The HDU list whose WCSCORR table should be appended to (the WCSCORR HDU
        must already exist)
    source : HDUList, optional
        The HDU list containing the extension from which to extract the WCS
        keywords to add to the WCSCORR table.  If None, the dest is also used
        as the source.
    extname : str, optional
        The extension name from which to take new WCS keywords.  If there are
        multiple extensions with that name, rows are added for each extension
        version.
    wcs_id : str, optional
        The name of the WCS to add, as in the WCSNAMEa keyword.  If
        unspecified, all the WCSs in the specified extensions are added.
    """

    if source is None:
        source = dest

    numext = fileutil.countExtn(source, extname)
    if numext == 0:
        raise ValueError('No %s extensions found in the source HDU list.'
                         % extname)

    # Current implementation assumes the same WCS keywords are in each
    # extension version; if this should not be assumed then this can be
    # modified...
    wcs_keys = altwcs.wcskeys(source[(extname, 1)].header)
    wcs_keys = filter(None, wcs_keys)
    wcshdr = stwcs.wcsutil.HSTWCS(source, ext=(extname, 1)).wcs2header()
    wcs_keywords = wcshdr.keys()

    if 'O' in wcs_keys:
        wcs_keys.remove('O') # 'O' is reserved for original OPUS WCS

    # If we're looking for a particular wcs_id, test ahead of time that it's
    # actually present in the specified extension headers
    if wcs_id:
        wcs_key = ''
        for wcs_key in wcs_keys:
            wcsname = source[(extname, 1)].header['WCSNAME' + wcs_key]
            if wcs_id == wcsname:
                break
        else:
            raise ValueError('A WCS with name %s was not found in the %s '
                             'extension headers in the source HDU list.'
                             % (wcs_id, extname))
        wcs_keys = [wcs_key] # We're only interested in this one

    # create new table for hdr and populate it with the newly updated values
    new_table = create_wcscorr(numrows=len(wcs_keys) * numext)
    old_table = dest['WCSCORR']

    idx = -1
    for wcs_key in wcs_keys:
        for extver in range(1, numext + 1):
            extn = (extname, extver)
            hdr = source[extn].header
            wcsname = hdr['WCSNAME' + wcs_key]
            selection = {'WCS_ID': wcsname, 'EXTVER': extver,
                         'WCS_key': wcs_key}

            # Ensure that an entry for this WCS is not already in the dest
            # table; if so just skip it
            rowind = find_wcscorr_row(old_table.data, selection)
            if np.any(rowind):
                continue

            idx += 1

            wcs = stwcs.wcsutil.HSTWCS(source, ext=extn, wcskey=wcs_key)
            wcshdr = wcs.wcs2header()

            # Update selection column values
            for key, val in selection.iteritems():
                new_table.data.field(key)[idx] = val

            for key in wcs_keywords:
                if key in new_table.data.names:
                    new_table.data.field(key)[idx] = wcshdr[key + wcs_key]

            prihdr = source[0].header
            for key in DEFAULT_PRI_KEYS:
                if key in new_table.data.names and prihdr.has_key(key):
                    new_table.data.field(key)[idx] = prihdr[key]

    # If idx was never incremented, no rows were added, so there's nothing else
    # to do...
    if idx < 0:
        return

    # Now, we need to merge this into the existing table
    rowind = find_wcscorr_row(old_table.data, {'wcs_id':''})
    old_nrows = np.where(rowind)[0][0]
    new_nrows = new_table.data.shape[0]

    # check to see if there is room for the new row
    if (old_nrows + new_nrows) > old_table.data.shape[0]:
        pad_rows = 5 * new_nrows
        # if not, create a new table with 'pad_rows' new empty rows
        upd_table = pyfits.new_table(old_table.columns,
                                     nrows=old_table.data.shape[0] + pad_rows)
    else:
        upd_table = old_table

    # Now, add
    for name in old_table.columns.names:
        upd_table.data.field(name)[old_nrows:old_nrows + new_nrows] = \
                new_table.data.field(name)
    upd_table.header.update('TROWS', old_nrows + new_nrows)

    # replace old extension with newly updated table extension
    dest['WCSCORR'] = upd_table


def restore_file_from_wcscorr(image, id='OPUS', wcskey=''):
    """ Copies the values of the WCS from the WCSCORR based on ID specified by user.
    The default will be to restore the original OPUS-derived values to the Primary WCS.
    If wcskey is specified, the WCS with that key will be updated instead.
    """

    if not isinstance(image, pyfits.HDUList):
        fimg = pyfits.open(image, mode='update')
        close_image = True
    else:
        fimg = image
        close_image = False
    numsci = fileutil.countExtn(fimg)
    wcs_table = fimg['WCSCORR']
    orig_rows = (wcs_table.data.field('WCS_ID') == 'OPUS')
    # create an HSTWCS object to figure out what WCS keywords need to be updated
    wcsobj = stwcs.wcsutil.HSTWCS(fimg,ext=('sci',1))
    wcshdr = wcsobj.wcs2header()
    for extn in range(1,numsci+1):
        # find corresponding row from table
        ext_rows = (wcs_table.data.field('EXTVER') == extn)
        erow = np.where(np.logical_and(ext_rows,orig_rows))[0][0]
        for key in wcshdr:
            if key in wcs_table.data.names: # insure that keyword is column in table
                tkey = key

                if 'orient' in key.lower():
                    key = 'ORIENT'
                if wcskey == '':
                    skey = key
                else:
                    skey = key[:7]+wcskey
                fimg['sci',extn].header.update(skey,wcs_table.data.field(tkey)[erow])
        for key in DEFAULT_PRI_KEYS:
            if key in wcs_table.data.names:
                if wcskey == '':
                    pkey = key
                else:
                    pkey = key[:7]+wcskey
                fimg[0].header.update(pkey,wcs_table.data.field(key)[erow])

    # close the image now that the update has been completed.
    if close_image:
        fimg.close()


def create_wcscorr(descrip=False, numrows=1, padding=0):
    """
    Return the basic definitions for a WCSCORR table.
    The dtype definitions for the string columns are set to the maximum allowed so
    that all new elements will have the same max size which will be automatically
    truncated to this limit upon updating (if needed).

    The table is initialized with rows corresponding to the OPUS solution
    for all the 'SCI' extensions.
    """

    trows = numrows + padding
    # define initialized arrays as placeholders for column data
    # TODO: I'm certain there's an easier way to do this... for example, simply
    # define the column names and formats, then create an empty array using
    # them as a dtype, then create the new table from that array.
    def_float64_zeros = np.array([0.0] * trows, dtype=np.float64)
    def_float64_ones = def_float64_zeros + 1.0
    def_float_col = {'format': 'D', 'array': def_float64_zeros.copy()}
    def_float1_col = {'format': 'D', 'array':def_float64_ones.copy()}
    def_str40_col = {'format': '40A',
                     'array': np.array([''] * trows, dtype='S40')}
    def_str24_col = {'format': '24A',
                     'array': np.array([''] * trows, dtype='S24')}
    def_int32_col = {'format': 'J',
                     'array': np.array([0]*trows,dtype=np.int32)}

    # If more columns are needed, simply add their definitions to this list
    col_names = [('CRVAL1',   def_float_col), ('CRVAL2',    def_float_col),
                 ('CRPIX1',   def_float_col), ('CRPIX2',    def_float_col),
                 ('CD1_1',    def_float_col), ('CD1_2',     def_float_col),
                 ('CD2_1',    def_float_col), ('CD2_2',     def_float_col),
                 ('CTYPE1',   def_str24_col), ('CTYPE2',    def_str24_col),
                 ('ORIENTAT', def_float_col), ('PA_V3',     def_float_col),
                 ('Delta_RA', def_float_col), ('Delta_Dec', def_float_col),
                 ('RMS_RA',   def_float_col), ('RMS_Dec',   def_float_col),
                 ('Delta_Orientat', def_float_col),
                 ('Delta_Scale', def_float1_col),
                 ('NMatch',   def_int32_col), ('Catalog',   def_str40_col)]

    # Define selector columns
    id_col = pyfits.Column(name='WCS_ID', format='40A',
                           array=np.array(['OPUS'] * numrows + [''] * padding,
                                          dtype='S24'))
    extver_col = pyfits.Column(name='EXTVER', format='I',
                               array=np.array(range(1, numrows + 1),
                                              dtype=np.int16))
    wcskey_col = pyfits.Column(name='WCS_key', format='A',
                               array=np.array(['O'] * numrows + [''] * padding,
                                              dtype='S'))
    # create list of remaining columns to be added to table
    col_list = [id_col, extver_col, wcskey_col] # start with selector columns
    for c in col_names:
        cdef = copy.deepcopy(c[1])
        col_list.append(pyfits.Column(name=c[0], format=cdef['format'],
                        array=cdef['array']))

    if descrip:
        col_list.append(
            pyfits.Column(name='Descrip', format='128A',
                          array=np.array(
                              ['Original WCS computed by OPUS'] * numrows,
                              dtype='S128')))

    # Now create the new table from the column definitions
    newtab = pyfits.new_table(pyfits.ColDefs(col_list), nrows=trows)
    # The fact that setting .name is necessary should be considered a bug in
    # pyfits.
    # TODO: Make sure this is fixed in pyfits, then remove this
    newtab.name = 'WCSCORR'

    return newtab

