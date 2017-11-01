# Driver functions to use AstrometryDB Restful service
# Based on interfaces/code provided by B. McLean 11-Oct-2017
import os
import time
import requests
from io import BytesIO
from lxml import etree

from stsci.tools import fileutil as fu
from stwcs.wcsutil import headerlet

import logging
logger = logging.getLogger('stwcs.updatewcs.astrometry_utils')


class AstrometryDB(object):

    serviceLocation='https://mastdev.stsci.edu/portal/astrometryDB/'
    headers = {'Content-Type': 'text/xml'}    

    available = True
    available_code = {'code' : "", 'text' : ""}


    def __init__(self, url=None, raise_errors=None):
        """
        Parameters
        ==========
        url : str
            User-provided URL for astrometry dB web-service to replaced
            default web-service
            
        raise_errors : bool, optional
             User can specify whether or not to turn off raising exceptions
             upon errors in getting or applying new astrometry solutions.  
             This will override the environment variable 'RAISE_PIPELINE_ERRORS'
             if set.
        """
        if url is not None:
            self.serviceLocation = url
        
        self.isAvailable() # determine whether service is available
        
        #
        # Implement control over behavior for error conditions
        # User provided input will always take precedent
        # Environment variable will also be recognized if no user-variable set
        # otherwise, it will turn off raising Exceptions
        #
        if raise_errors is not None:
            self.raise_errors = raise_errors 
        elif 'RAISE_PIPELINE_ERRORS' in os.environ:
            self.raise_errors = os.environ['RAISE_PIPELINE_ERRORS']
        else:
            self.raise_errors = False

    def updateObs(self, obsname):
        """ Method to update observation with any available solutions
        
        Parameters
        ==========
        obsname : str
           Filename (without path?) for observation to be updated            

        """
        if not self.available:
            logger.warning("AstrometryDB not available. NO Updates performed for {}".format(obsname))
            return
        
        obspath, obsroot = os.path.split(obsname)
        # use `obspath` for location of output files, if anything gets written out
        observationID = fu.getRootname(obsroot)
        
        headerlets = self.getObservation(observationID)
        #
        # apply to file...
        fileobj = pf.open(obsname, mode='update')
        for h in headerlets.keys():
            h.apply_as_alternate(fileobj, wcsname=h, attach=True)    
               
        fileobj.close()
        
    def getObservation(self,observationID):
        """ Method for getting solutions for observation from AstrometryDB
        
        Parameters
        ==========
        observationID : str
            base rootname for observation to be updated (eg., `iab001a1q`)
            
        """
        if not self.available:
            logger.warning("AstrometryDB not available. NO Updates performed for {}".format(observationID))
            return None

        serviceEndPoint=self.serviceLocation+'observation/read/'+observationID

        try:
            logger.info('Accessing AstrometryDB service : {}'.format(serviceEndPoint))
            r = requests.get(serviceEndPoint,headers=self.headers)
            if r.status_code == requests.codes.ok:
                logger.info('AstrometryDB service call succeeded') 
            else:
                logger.warning(" AstrometryDB service call failed")
                logger.warning("    Status: {}".format(r.status_code))
                logger.warning("    Error:  {}".format(r.text))
                return None
        except:
            logger.warning('AstrometryDB service call failed')
            logger.warning("    Status: {}".format(r.status_code))
            return None
        
        # Now, interpret return value for observation into separate headerlets
        # to be appended to observation
        headerlets = {}
        tree = BytesIO(r.content)
        solutions = []
        # get names of solutions in database
        for _,element in etree.iterparse(tree, tag='solution'):
            solutions.append(element[1].text)
        # Now use these names to get the actual updated solutions
        headers = {'Content-Type': 'application/fits'}
        for solutionID in solutions:
            serviceEndPoint=serviceLocation+'observation/read/'+observationID+'?wcsname='+solutionID
            r_solution = requests.get(serviceEndPoint,headers=headers)
            if r_solution.status_code == requests.codes.ok:
                headerlets[solutionID] = headerlet.Headerlet(BytesIO(r_solution.content).getvalue())
        
        return headerlets
        
    def isAvailable(self):
        """ Method to test availability of astrometryDB web-service
        """            
        serviceEndPoint=self.serviceLocation+'availability'

        try:
            r = requests.get(serviceEndPoint,headers=self.headers)
            if r.status_code == requests.codes.ok:
                logger.info('AstrometryDB service available...')
                self.avaiable_code['code'] = r.status_code
                self.available_code['text'] = 'Available'
            else:
                logger.warning('WARNING : AstrometryDB service unavailable!')
                logger.warning('    AstrometryDB status: {}'.format(r.status_code))
                logger.warning('    AstrometryDB text: {}'.format(r.text))
                self.available_code['code'] = r.status_code
                self.available_code['text'] = r.text
                self.available=False
        except:
            logger.warning('WARNING : AstrometryDB service inaccessible!')
            logger.warning('    AstrometryDB status: {}'.format(r.status_code))
            self.available_code['code'] = r.status_code
            self.available=False
        

def apply_astrometric_updates(obsnames, **pars):
    """ Functional stand-alone interface for applying new astrometric solutions
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
            retrieving a solution from the database or applying the new solution
            to the observation. If None, it will look to see whether the 
            environmental variable `RAISE_PIPELINE_ERRORS` was set, otherwise,
            it will default to 'False'. 
    """
    
    if type(obsnames) != type([]):
        obsnames = [obsnames]
        
    url = pars.get('url',None)
    raise_errors = pars.get('raise_errors', None)
    
    db = AstrometryDB(url=url, raise_errors=raise_errors)
    for obs in obsnames:
        db.updateObs(obs)
        

    
