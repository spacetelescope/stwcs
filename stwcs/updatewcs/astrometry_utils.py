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

ASTROMETRY_STEP_CONTROL - String specifying whether or not to perform the
                          astrometry update processing at all.
                          Valid Values: "ON", "On", "on", "OFF", "Off", "off"
                          If not set, default value is "ON".
"""
import os
import atexit

import requests
from io import BytesIO
from lxml import etree

from astropy.io import fits as pf

from stwcs.wcsutil import headerlet

import logging
logger = logging.getLogger('stwcs.updatewcs.astrometry_utils')

atexit.register(logging.shutdown)

# Definitions of environment variables used by this step
astrometry_db_envvar = "ASTROMETRY_SERVICE_URL"
pipeline_error_envvar = "RAISE_PIPELINE_ERRORS"
astrometry_control_envvar = "ASTROMETRY_STEP_CONTROL"


class AstrometryDB(object):
    """Base class for astrometry database interface."""

    serviceLocation = 'https://mastdev.stsci.edu/portal/astrometryDB/'
    headers = {'Content-Type': 'text/xml'}

    available = True
    available_code = {'code': "", 'text': ""}

    def __init__(self, url=None, raise_errors=None, perform_step=True,
                 write_log=False):
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

        perform_step : bool, optional
            Specify whether or not to perform this step.  This will
            completely override the setting of the `ASTROMETRY_STEP_CONTROL`
            environment variable.  Default value: True.

        write_log : bool, optional
            Specify whether or not to write a log file during processing.
            Default: False

        """
        self.perform_step = perform_step
        # Check to see whether an environment variable has been set
        if astrometry_control_envvar in os.environ:
            val = os.environ[astrometry_control_envvar].lower()
            if val == 'off':
                self.perform_step = False
            elif val == 'on':
                self.perform_step = True
            else:
                l = "Environment variable {} not set to valid value".\
                    format(astrometry_control_envvar)
                l += "\t Valid values: on or off (case-insensitive)"
                raise ValueError(l)
            logger.info("Astrometry step operation set to {}".
                        format(self.perform_step))
        if not self.perform_step:
            logger.info("Astrometry update step has been turned off")
            logger.info("\tNo updates will be performed!")
            return

        if write_log:
            formatter = logging.Formatter(
                        "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            log_filename = 'astrometry.log'
            fh = logging.FileHandler(log_filename, mode='w')
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(formatter)
            logger.addHandler(fh)
        logger.setLevel(logging.INFO)

        # check to see whether any URL has been specified as an
        # environmental variable.
        if astrometry_db_envvar in os.environ:
            self.serviceLocation = os.environ[astrometry_db_envvar]

        if url is not None:
            self.serviceLocation = url
        #
        # Implement control over behavior for error conditions
        # User provided input will always take precedent
        # Environment variable will also be recognized if no user-variable set
        # otherwise, it will turn off raising Exceptions
        #
        self.raise_errors = False
        if raise_errors is not None:
            self.raise_errors = raise_errors
            logger.info("Setting `raise_errors` to {}".format(raise_errors))
        if pipeline_error_envvar in os.environ:
            self.raise_errors = True

        self.isAvailable()  # determine whether service is available

        # Initialize attribute to keep track of type of observation
        self.new_observation = False

    def updateObs(self, obsname):
        """Update observation with any available solutions.

        Parameters
        ==========
        obsname : str
           Filename for observation to be updated

        """
        if not self.perform_step:
            return

        # Parse observation name
        obspath, obsroot = os.path.split(obsname)
        # use `obspath` for location of output files,
        #    if anything gets written out
        observationID = obsroot.split('_')[:1][0]
        logger.info("Updating astrometry for {}".format(observationID))
        #
        # apply to file...
        fileobj = pf.open(obsname, mode='update')

        # take inventory of what hdrlets are already appended to this file
        hdrnames = headerlet.get_headerlet_kw_names(fileobj, 'hdrname')

        headerlets, best_solution_id = self.getObservation(observationID)
        if headerlets is None:
            logger.warning("Problems getting solutions from database")
            logger.warning(" NO Updates performed for {}".format(
                           observationID))
            if self.raise_errors:
                raise ValueError("No new solution found in AstrometryDB.")
            else:
                return

        # If no headerlet found in database, update database with this WCS
        if self.new_observation:
            logger.warning(" No new solution found in AstrometryDB.")
            logger.warning(" Updating database with initial WCS {}".
                           format(observationID))
            hlet_buffer = BytesIO()
            hlet_new = headerlet.create_headerlet(fileobj)
            newhdrname = hlet_new[0].header['hdrname']
            hlet_new.writeto(hlet_buffer)

            logger.info("Updating AstrometryDB with entry for {}".format(
                        observationID))
            logger.info("\t using WCS with HDRNAME={}".format(newhdrname))
            # Add WCS solution from this observation to the database
            self.addObservation(observationID, hlet_buffer)
        else:
            # Attach new unique hdrlets to file...
            logger.info("Updating {} with:".format(observationID))
            for h in headerlets:
                newhdrname = headerlets[h][0].header['hdrname']
                if newhdrname in hdrnames:
                    continue  # do not add duplicate hdrlets
                # Add solution as an alternate WCS
                try:
                    if best_solution_id and newhdrname == best_solution_id:
                        # replace primary WCS with this solution
                        headerlets[h].apply_as_primary(fileobj)
                        logger.info('Replacing primary WCS with')
                        logger.info('\tHeaderlet with HDRNAME={}'.format(
                                     newhdrname))
                    else:
                        logger.info("\tHeaderlet with HDRNAME={}".format(
                                    newhdrname))
                        headerlets[h].attach_to_file(fileobj)
                except ValueError:
                    pass

        fileobj.close()

    def findObservation(self, observationID):
        """Find whether there are any entries in the AstrometryDB for
        the observation with `observationID`.

        Parameters
        ==========
        observationID : str
            base rootname for observation to be updated (eg., `iab001a1q`)

        Return
        ======
        entry : obj
            Database entry for this observation, if found.
            It will return None if there was an error in accessing the
            database and `self.raise_errors` was not set to True.
        """
        if not self.perform_step:
            return None

        serviceEndPoint = self.serviceLocation + \
            'observation/read/' + observationID

        try:
            logger.info('Accessing AstrometryDB service :')
            logger.info('\t{}'.format(serviceEndPoint))
            r = requests.get(serviceEndPoint, headers=self.headers)
            if r.status_code == requests.codes.ok:
                logger.info('AstrometryDB service call succeeded')
            elif r.status_code == 404:
                # This code gets returned if exposure is not found in database
                # Never fail for this case since all new observations
                # will result in this error
                logger.info("No solutions found in database for {}".
                            format(observationID))
                self.new_observation = True
            else:
                logger.warning(" AstrometryDB service call failed")
                logger.warning("    Status: {}".format(r.status_code))
                logger.warning("    {}".format(r.reason))
                if self.raise_errors:
                    e = "AstrometryDB service could not be connected!"
                    raise requests.RequestException(e)
                else:
                    return None
        except Exception:
            logger.warning('AstrometryDB service call failed')
            logger.warning("    Status: {}".format(r.status_code))
            logger.warning("    {}".format(r.reason))

            if self.raise_errors:
                l = 'AstrometryDB service call failed with reason:\n\t"{}"'.\
                    format(r.reason)
                l += '\n\tStatus code = {}'.format(r.status_code)
                raise requests.RequestException(l)
            else:
                return None
        return r

    def getObservation(self, observationID):
        """Get solutions for observation from AstrometryDB.

        Parameters
        ==========
        observationID : str
            base rootname for observation to be updated (eg., `iab001a1q`)

        Return
        ======
        headerlets : dict
            Dictionary containing all solutions found for exposure in the
            form of headerlets labelled by the name given to the solution in
            the database.
        """
        if not self.perform_step:
            return None, None

        if not self.available:
            logger.warning("AstrometryDB not available.")
            logger.warning("NO Updates performed for {}".format(observationID))
            if self.raise_errors:
                raise ConnectionError("AstrometryDB not accessible.")
            else:
                return None, None

        r = self.findObservation(observationID)

        if r is None or self.new_observation:
            return r, None
        else:
            # Now, interpret return value for observation into separate
            # headerlets to be appended to observation
            headerlets = {}
            tree = BytesIO(r.content)
            solutions = []
            # get names of solutions in database
            for _, element in etree.iterparse(tree, tag='solution'):
                s = element[1].text
                if s:
                    solutions.append(s)
            # get name of best solution specified by database
            tree.seek(0)
            best_solution_id = None
            for _, element in etree.iterparse(tree, tag='bestsolutionid'):
                best_solution_id = element
                break
            if best_solution_id == '':
                best_solution_id = None

            # Now use these names to get the actual updated solutions
            headers = {'Content-Type': 'application/fits'}
            for solutionID in solutions:
                serviceEndPoint = self.serviceLocation + \
                    'observation/read/' + observationID + \
                    '?wcsname='+solutionID
                r_solution = requests.get(serviceEndPoint, headers=headers)
                if r_solution.status_code == requests.codes.ok:
                    hlet_bytes = BytesIO(r_solution.content).getvalue()
                    hlet = headerlet.Headerlet(file=hlet_bytes)
                    hlet.init_attrs()
                    if hlet[0].header['hdrname'] == 'OPUS':
                        hdrdate = hlet[0].header['date'].split('T')[0]
                        hlet[0].header['hdrname'] += hdrdate
                    headerlets[solutionID] = hlet

            if not solutions:
                logger.warning("No new WCS's found for {}".format(observationID))
                logger.warning("No updates performed...")
                
            return headerlets, best_solution_id

    def addObservation(self, observationID, new_solution):
        """Add WCS from current observation to database"""
        if not self.perform_step:
            return

        serviceEndPoint = self.serviceLocation+'observation/create'
        headers = {'Content-Type': 'application/octet-stream'}

        r = requests.post(serviceEndPoint, data=new_solution, headers=headers)
        if r.status_code == requests.codes.ok:
            logger.info("AstrometryDB service updated with new entry for {}".
                        format(observationID))
        else:
            l = "Problem encountered when adding {} to database".\
                           format(observationID)
            if self.raise_errors:
                raise Exception(l)
            else:
                logger.warning(l)

    def isAvailable(self):
        """Test availability of astrometryDB web-service."""
        if not self.perform_step:
            return

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
                if self.raise_errors:
                    e = "AstrometryDB service unavailable!"
                    raise ConnectionRefusedError(e)

        except Exception as err:
            logger.warning('WARNING : AstrometryDB service inaccessible!')
            logger.warning('    AstrometryDB called: {}'.format(
                                self.serviceLocation))
            self.available = False
            if self.raise_errors:
                raise ConnectionError from err


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
    if not isinstance(obsnames, list):
        obsnames = [obsnames]

    url = pars.get('url', None)
    raise_errors = pars.get('raise_errors', None)

    db = AstrometryDB(url=url, raise_errors=raise_errors)
    for obs in obsnames:
        db.updateObs(obs)
