
Introduction
============

Calibration of the HST Advanced Camera for Surveys (HST/ACS) distortion requires the use 
of several components to the distortion correction; namely, polynomial coefficients, a correction 
for velocity aberration, a time-dependent skew, a lookup table for non-polynomial terms, and a detector defect 
correction. Each of these terms has been derived as part of the calibration effort to address 
separate aspects of the distortion that affects ACS observations. Ideally, each would be applied
independently in the same manner used for deriving the original calibration reference information, 
with the time-dependent skew being folded into the other terms. However, the software for 
applying the distortion models does not support this option. In fact, there is no clear 
accepted standard for specifying distortion corrections in FITS headers. Instead, there are 
several separate proposals for specifying aspects of the distortion, but none by themselves 
allows us to fully specify the distortion already calibrated for ACS, let alone in a modular, 
efficient manner.

This paper describes a composite implementation of a select set of proposed standards which 
supports all aspects of the distortion models for HST instruments without breaking any of the 
conventions/standards. The rules for merging the proposed standards allow software to be defined 
to apply each aspect of this proposal as a separate option while defining the requirements 
necessary to allow them to work together when specified in the header. As a result, the separate 
components essentially become tools where only those conventions appropriate to the observation 
can be used as needed. 


