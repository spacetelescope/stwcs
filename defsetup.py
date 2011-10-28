from __future__ import division # confidence high

import sys

pkg =  ["stwcs", 'stwcs.updatewcs', 'stwcs.wcsutil', 'stwcs.distortion', 'stwcs.gui']
pkg_gui = "stwcs/gui"

setupargs = {
    'version' :         "0.8",
    'description' :		"Recomputes the WCS of an HST observation and puts all istortion corrections in the headers.",
    'package_dir': {'stwcs':'lib/stwcs', 'stwcs.updatewcs': 'lib/stwcs/updatewcs', 'stwcs.gui':'lib/stwcs/gui', 
                    'stwcs.wcsutil': 'lib/stwcs/wcsutil', 'stwcs.distortion': 'lib/stwcs/distortion'},
    'data_files' :        [( pkg_gui+"/pars", ['lib/stwcs/gui/pars/*']),
                            ( pkg_gui+"/htmlhelp", ['lib/stwcs/gui/htmlhelp/*.html']),
                            ( pkg_gui, ['lib/stwcs/gui/*.help'])],
    'author' :		    "Nadia Dencheva, Warren Hack",
    'author_email' :    "help@stsci.edu",
    'license' :		    "http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE",
    'platforms' :	    ["Linux","Solaris","Mac OS X", "Windows"],
}

