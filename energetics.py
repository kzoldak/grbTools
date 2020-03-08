from __future__ import division # just incase this is running on Python 2.7
from collections import OrderedDict
from scipy import integrate
from math import pi
import numpy as np

from numpy import exp, log, log10
from math import atanh

#from .models import *

from grbTools.cosmology.luminositydistance import LumDist  # basic version



__all__ = ["lpow", "copl", "band", "grbm", "sbpl", "bbody", 
           "calc_flux", "calc_fluence", "calc_eiso"]



class Energetics(object):
    def __init__(self, modelname=None, parameters=None, duration=None, 
        emin=None, emax=None, redshift=None, program=None):
        if program == 'general'
        pass



    def calc_flux(self):
        pass

    def calc_fluence(self, modelname, parameters, duration, emin, emax, redshift=None, 
        photonflux=False, program='general'):
        flux = calc_flux(modelname, parameters, emin, emax, redshift, photonflux)
        fluence = flux * duration
        return fluence 

    def calc_eiso(modelname, parameters, duration, emin, emax, redshift, program='general'):
        '''
        Notes:
        ------
        Luminosity Distance Calculation can be found here:
        from grbTools.Cosmology.luminositydistance import LumDist
        Import LumDist and read its documentation to learn more on the 
        cosmology we use. 

        '''
        Energetics.__init__(modelname, parameters, duration, emin, emax, redshift, program='general')
        
        fluence = self.calc_fluence(modelname, parameters, duration, emin, emax, 
            redshift, photonflux=False)
        DL = LumDist(redshift)
        eiso = (4.0*pi*pow(DL,2)*fluence)/(1.+redshift)
        return eiso

    def calc_flux_uncertainty(self, modelname, program):
        if program == 'general':
            raise Exception('No propagation of uncertainty for general models.')
        pass

    def calc_fluence_uncertainty(self, modelname, program):
        if program == 'general':
            raise Exception('No propagation of uncertainty for general models.')
        pass

    def calc_eiso_uncertainty(self, modelname, program):
        if program == 'general':
            raise Exception('No propagation of uncertainty for general models.')
        pass

    def calc_epeak_uncertainty(self, modelname, program):
        if program == 'general':
            raise Exception('No propagation of uncertainty for general models.')
        pass



