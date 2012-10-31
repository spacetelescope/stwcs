FITS Distortion Proposal
=========================

The current FITS Distortion Paper conventions [DistortionPaper]_ provide a mechanism for specifying either a lookup table 
or polynomial model for the distortion of each axis. The standard states in Section 2.1:

 Note that the prior distortion functions, :math:`\delta_p(p)`, operate on pixel coordinates (i.e. 
 :math:`p` rather than :math:`p-r`), and that the independent variables of the distortion functions 
 are the *uncorrected* pixel or intermediate pixel coordinates. That is, for example, 
 we do not allow the possibility of

.. math::
   :label: Equation 1

   q'_{3} = q_{3} + \delta_{q_{3}}(q'_{1},q'_{2})

The keywords used for describing these corrections use the syntax given in Table 2 of the FITS Distortion Paper. 
For our purposes, the keywords of interest are those related to lookup tables; namely, 

::

 CPDISja        string    2.4.1 distortion code new Prior distortion function type.
 DPja           record    2.4.2 distortion parameter new Parameter for a prior distortion 
                                  function, for use in an image header
                          
This syntax only provides the option to specify one correction at a time for each 
axis of the image. This precludes being able to use this convention to specify both 
a lookup table and a polynomial model at the same time for the same axis. It does not 
state what should be done if the polynomial has been specified using a different 
convention, for example, the SIP convention. Thus, SIP and FITS Distortion Paper should not be 
seen as mutually exclusive. In fact, they may work together rather naturally since the 
SIP and FITS Distortion Paper conventions both assume the corrections will work on the input pixel 
and add to the output frame. 

The sample header in :ref:`appendix1` shows how these keywords get populated for
an actual reference file; specifically, an NPOLFILE as described in the next section.


