import sys

pkg =  ["updatewcs", 'updatewcs.hstwcs', 'updatewcs.wcsutil', 'updatewcs.distortion']

setupargs = {
    'version' :         "0.1",
    'description' :		"Recomputes the WCS of an HST observation and puts all istortion corrections in the headers.",
    'package_dir': {'updatewcs':'lib', 'updatewcs.hstwcs': 'hstwcs',
                    'updatewcs.wcsutil': 'wcsutil', 'updatewcs.distortion': 'distortion'},

    'author' :		    "Nadia Dencheva, Warren Hack",
    'author_email' :    "help@stsci.edu",
    'license' :		    "http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE",
    'platforms' :	    ["Linux","Solaris","Mac OS X", "Windows"],
}

