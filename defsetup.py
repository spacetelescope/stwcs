import sys

pkg =  ["hstwcs", 'hstwcs.updatewcs', 'hstwcs.wcsutil', 'hstwcs.distortion']

setupargs = {
    'version' :         "0.1",
    'description' :		"Recomputes the WCS of an HST observation and puts all istortion corrections in the headers.",
    'package_dir': {'hstwcs':'lib', 'hstwcs.updatewcs': 'updatewcs',
                    'hstwcs.wcsutil': 'wcsutil', 'hstwcs.distortion': 'distortion'},

    'author' :		    "Nadia Dencheva, Warren Hack",
    'author_email' :    "help@stsci.edu",
    'license' :		    "http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE",
    'platforms' :	    ["Linux","Solaris","Mac OS X", "Windows"],
}

