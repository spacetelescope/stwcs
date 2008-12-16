

# This dictionary maps an instrument into an instrument class
# The instrument class handles instrument specific keywords

inst_mappings={'WFPC2': 'WFPC2WCS',
                'ACS': 'ACSWCS'
                }


# A list of instrument specific keywords
# Every instrument class must have methods which define each of these
# as class attributes.
ins_spec_kw = ['ltv1', 'ltv2', 'parity', 'binned','vafactor', 'chip', 
            'naxis1', 'naxis2', 'filter1', 'filter2', 'detector']

# A list of keywords defined in the primary header.
# The HSTWCS class sets this as attributes 
prim_hdr_kw = ['detector', 'offtab', 'idctab', 'date-obs', 
              'pa_v3', 'ra_targ', 'dec_targ']

# These are the keywords which are archived before MakeWCS is run
basic_wcs = ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CRVAL1','CRVAL2','CTYPE1', 'CTYPE2',
            'CRPIX1', 'CRPIX2', 'CTYPE1', 'CTYPE2', 'ORIENTAT', 'NAXIS1', 'NAXIS2']
            
