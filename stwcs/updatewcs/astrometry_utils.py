"""Interface for Astrometry database service.

This module contains interface functions for the AstrometryDB Restful service
based on interfaces/code provided by B. McLean 11-Oct-2017.

The code checks for the existence of environmental variables to help control
the operation of this interface; namely,

RAISE_PIPELINE_ERRORS - boolean to specify whether to raise exceptions
                        during processing or simply log errors and quit
                        gracefully.  If not set, default behavior will be
                        to log errors and quit gracefully.

ASTROMETRY_SERVICE_URL - URL pointing to user-specified web service that
                        will provide updated astrometry solutions for
                        observations being processed.  This will replace
                        the built-in URL included in the base class.  This
                        value will also be replaced by any URL provided by
                        the user as an input parameter `url`.

"""
import os
import requests
from io import BytesIO
from lxml import etree

from astropy.io import fits as pf

from stsci.tools import fileutil as fu
from stwcs.wcsutil import headerlet

import logging
logger = logging.getLogger('stwcs.updatewcs.astrometry_utils')

astrometry_db_envvar = "ASTROMETRY_SERVICE_URL"
pipeline_error_envvar = "RAISE_PIPELINE_ERRORS"


class AstrometryDB(object):
    """Base class for astrometry database interface."""

    serviceLocation = 'https://mastdev.stsci.edu/portal/astrometryDB/'
    headers = {'Content-Type': 'text/xml'}

    available = True
    available_code = {'code': "", 'text': ""}

    def __init__(self, url=None, raise_errors=None):
        """Initialize class with user-provided URL.

        Parameters
        ==========
        url : str
            User-provided URL for astrometry dB web-service to replaced
            default web-service.  Any URL specified here will override
            any URL specified as the environment variable
            `ASTROMETRY_SERVICE_URL`.  This parameter value, if specified,
            and the environmental variable will replace the built-in default
            URL included in this class.

        raise_errors : bool, optional
             User can specify whether or not to turn off raising exceptions
             upon errors in getting or applying new astrometry solutions.
             This will override the environment variable
             'RAISE_PIPELINE_ERRORS' if set.

        """
        # check to see whether any URL has been specified as an
        # environmental variable.
        if astrometry_db_envvar in os.environ:
            self.serviceLocation = os.environ[astrometry_db_envvar]

        if url is not None:
            self.serviceLocation = url

        self.isAvailable()  # determine whether service is available

        #
        # Implement control over behavior for error conditions
        # User provided input will always take precedent
        # Environment variable will also be recognized if no user-variable set
        # otherwise, it will turn off raising Exceptions
        #
        if raise_errors is not None:
            self.raise_errors = raise_errors
        elif pipeline_error_envvar in os.environ:
            self.raise_errors = os.environ[pipeline_error_envvar]
        else:
            self.raise_errors = False

    def updateObs(self, obsname):
        """Update observation with any available solutions.

        Parameters
        ==========
        obsname : str
           Filename (without path?) for observation to be updated

        """
        if not self.available:
            logger.warning("AstrometryDB not available.")
            logger.warning(" NO Updates performed for {}".format(obsname))
            return

        obspath, obsroot = os.path.split(obsname)
        # use `obspath` for location of output files,
        #    if anything gets written out
        observationID = fu.getRootname(obsroot)

        headerlets = self.getObservation(observationID)
        #
        # apply to file...
        fileobj = pf.open(obsname, mode='update')
        for h in headerlets.keys():
            h.apply_as_alternate(fileobj, wcsname=h, attach=True)

        fileobj.close()

    def getObservation(self, observationID):
        """Get solutions for observation from AstrometryDB.

        Parameters
        ==========
        observationID : str
            base rootname for observation to be updated (eg., `iab001a1q`)

        """
        if not self.available:
            logger.warning("AstrometryDB not available.")
            logger.warning(" NO Updates performed for {}".format(
                           observationID))
            return None

        serviceEndPoint = self.serviceLocation + \
            'observation/read/' + observationID

        try:
            logger.info('Accessing AstrometryDB service : {}'.format(
                        serviceEndPoint))
            r = requests.get(serviceEndPoint, headers=self.headers)
            if r.status_code == requests.codes.ok:
                logger.info('AstrometryDB service call succeeded')
            else:
                logger.warning(" AstrometryDB service call failed")
                logger.warning("    Status: {}".format(r.status_code))
                logger.warning("    Error:  {}".format(r.text))
                return None
        except Exception:
            logger.warning('AstrometryDB service call failed')
            logger.warning("    Status: {}".format(r.status_code))
            return None

        # Now, interpret return value for observation into separate headerlets
        # to be appended to observation
        headerlets = {}
        tree = BytesIO(r.content)
        solutions = []
        # get names of solutions in database
        for _, element in etree.iterparse(tree, tag='solution'):
            solutions.append(element[1].text)
        # Now use these names to get the actual updated solutions
        headers = {'Content-Type': 'application/fits'}
        for solutionID in solutions:
            serviceEndPoint = self.serviceLocation + \
                'observation/read/' + observationID + '?wcsname='+solutionID
            r_solution = requests.get(serviceEndPoint, headers=headers)
            if r_solution.status_code == requests.codes.ok:
                headerlets[solutionID] = headerlet.Headerlet(
                    BytesIO(r_solution.content).getvalue())

        return headerlets

    def isAvailable(self):
        """Test availability of astrometryDB web-service."""
        serviceEndPoint = self.serviceLocation+'availability'

        try:
            r = requests.get(serviceEndPoint, headers=self.headers)

            if r.status_code == requests.codes.ok:
                logger.info('AstrometryDB service available...')
                self.available_code['code'] = r.status_code
                self.available_code['text'] = 'Available'
                self.available = True
            else:
                logger.warning('WARNING : AstrometryDB service unavailable!')
                logger.warning('          AstrometryDB called: {}'.format(
                                self.serviceLocation))
                logger.warning('          AstrometryDB status: {}'.format(
                               r.status_code))
                logger.warning('          AstrometryDB text: {}'.format(
                               r.text))
                self.available_code['code'] = r.status_code
                self.available_code['text'] = r.text
                self.available = False
        except Exception as err:
            logger.warning('WARNING : AstrometryDB service inaccessible!')
            logger.warning('    AstrometryDB called: {}'.format(
                                self.serviceLocation))
            logger.warning('    AstrometryDB status: {}'.format(r.status_code))
            self.available_code['code'] = r.status_code
            self.available = False


def apply_astrometric_updates(obsnames, **pars):
    """Apply new astrometric solutions to observation.

    Functional stand-alone interface for applying new astrometric solutions
    found in the astrometry dB to the given observation(s).

    Parameters
    ==========
    obsnames : str, list
        Filename or list of filenames of observation(s) to be updated

    url : str, optional
        URL of astrometry database web-interface to use.
        If None (default), it will use built-in URL for STScI web-interface

    raise_errors : bool, optional
        Specify whether or not to raise Exceptions and stop processing
        when an error in either accessing the database,
        retrieving a solution from the database or applying the new
        solution to the observation. If None, it will look to see whether
        the environmental variable `RAISE_PIPELINE_ERRORS` was set,
        otherwise, it will default to 'False'.

    """
    if isinstance(obsnames, type([])):
        obsnames = [obsnames]

    url = pars.get('url', None)
    raise_errors = pars.get('raise_errors', None)

    db = AstrometryDB(url=url, raise_errors=raise_errors)
    for obs in obsnames:
        db.updateObs(obs)
