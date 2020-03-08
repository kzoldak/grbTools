from __future__ import division
from collections import OrderedDict
from scipy import integrate
from math import pi
import numpy as np

from .models import *

class XspecEnergetics(object):
    """
    Parameters:
    -----------
    model : str, model name. Full name with + between additive models.
    parameters : ordered dict, see below.
    emin : float, min energy integrand.
    emax : float, max energy integrand.
    redshift : float, redshift.
    photonflux : True or False, default is False. True returns photon fluxes 
                instead of energy fluxes. 
                    
    This class is only to be used when reading directly from XSPEC files 
    returned by the program after spectral fitting. Therefore, we excect 
    the user to know what order the parameters should be in and how to 
    properly construct the parameter dictionaries needed to calculate 
    fluxes. Therefore we don't provide a function for returning an example 
    of one for each model. We provide a few examples here.

    Here are a few examples:
    ------------------------

    from collections import OrderedDict
    
    # For the GRBM model.
    pars    = OrderedDict()
    pars['grbm']     = OrderedDict()
    pars['grbm']['alpha__1']= -1.0
    pars['grbm']['beta__2']= -2.3
    pars['grbm']['tem__3']= 300.0
    pars['grbm']['norm__4']= 0.01

    # For the GRBM+LPOW model.
    pars    = OrderedDict()
    pars['grbm']     = OrderedDict()
    pars['grbm']['alpha__1']= -1.0
    pars['grbm']['beta__2']= -2.3
    pars['grbm']['tem__3']= 300.0
    pars['grbm']['norm__4']= 0.01
    pars['lpow']     = OrderedDict()
    pars['lpow']['plIndex__5']= -1.0
    pars['lpow']['norm__6']= 0.001

    Notes on these dictionaries:
    ----------------------------
    - Always provide nested dictionaries, one for each submodel.
    - Make sure you use an OrderedDict from collections. 
    - Structure the OrderedDict for each submodel in the order it was 
      defined in XSPEC. I.e., if your model is "grbm+bbody+lpow", 
      then the first model in the dict should be "grbm", the second "bbody", 
      and the last should be "lpow". The parameters for each submodel must 
      be placed in order (as the function takes them) as well. 
    - Although the first parameter in each model function is 'energy', 
      it is not one that should be included in parameter dicts and should 
      not be passed to the function. The integration takes care of that.

    Functions and their parameters are:
    -----------------------------------
        grbm(energy, alpha, beta, tem, norm)
        sbpl(energy, alpha, beta, ebreak, norm)
        cutoffpl(energy, PhoIndex, HighECut, norm)
        bbody(energy, kT, norm)
        lpow(energy, plIndex, norm)
        powerlaw(energy, PhoIndex, norm)

    """

    def __init__(self, modelname, parameters, emin, emax, redshift=None, photonflux=False):
        self.modelname = modelname
        self.parameters = parameters
        self.emin = emin
        self.emax = emax
        self.redshift = redshift
        self.photonflux = photonflux 
        self.model_options = 'grbm sbpl cutoffpl bbody lpow powerlaw'.split()
        
        if modelname not in self.model_options:
            msg = "'modelname' must be one of following: %r"%self.model_options
            raise Exception(msg)
          

    def calc_flux(self):
        """

        Parameters:
        -----------
        model : str, model name. Full name with + between additive models.
        parameters : ordered dict, see below.
        emin : float, min energy integrand.
        emax : float, max energy integrand.
        redshift : float, redshift.
        photonflux : True or False, default is False. True returns photon fluxes 
                    instead of energy fluxes. 

        
        This function calculates energy flux from an xspec model fit. 
        We build a dictionary of the model and parameters and then pass it 
        to Calc_Flux. This particular Calc_Flux function was set up for XSPEC 
        models only. 

        pars dictionary
        ---------------
        pars is an ordered dictionary, meaning the keys and values stay in the 
          order they were entered. This is critical for this program because the 
          Calc_Flux function was set up so that a list of parameters is passed 
          to the model functions, not a dictionary. 

          For example:
          grbm(energy, alpha, beta, tem, norm) 
          is the model set up for the band function. 

          lpow(plIndex, norm) 
          is the model set up for the power-law function. 

          Within the pars dictionary, the parameters must be set up as so:
              pars = OrderedDict()

              pars['grbm'] = OrderedDict()
              pars['grbm']['alpha__1']= -1.0
              pars['grbm']['beta__2']= -2.3
              pars['grbm']['tem__3']= 300.0
              pars['grbm']['norm__4']= 0.01

              pars['lpow'] = OrderedDict()
              pars['lpow']['plIndex__5']= -1.0
              pars['lpow']['norm__6']= 0.001

          Notice the order of the parameters. 
          If we do 
              pars['grbm']['beta__2']= -2.3
              pars['grbm']['alpha__1']= -1.0
              pars['grbm']['tem__3']= 300.0
              pars['grbm']['norm__4']= 0.01
          instead, then the band function reads the -2.3 as alpha, not beta. 
          Thus, having the parameter names defined within the dictionary is 
          not there to help the model functions (being integrated) know the 
          proper parameter values, but to keep things organized on our end. 
          
          You could do the following, and it wouldn't matter:
              pars['grbm']['alp']= -2.3
              pars['grbm']['bet']= -1.0
              pars['grbm']['ech']= 300.0
              pars['grbm']['amp']= 0.01
          But we don't recommend doing this. 


          Setting up models
          ------------------
          Regardless of whether or not you have additive models, you should 
          always define the model name as the first level of keys in the nested 
          dictionary. 
          For example, if we only have the band function results, we would still 
          do:
                pars['grbm'] = OrderedDict()
                pars['grbm']['alpha__1'] = -1.0
                pars['grbm']['beta__2'] = -2.3
                pars['grbm']['tem__3'] = 300.0
                pars['grbm']['norm__4'] = 0.01


        """
        keVtoerg = 1.60217657E-9

        if redshift is not None:
            emin = emin/(1.+redshift)
            emax = emax/(1.+redshift)
        else:
            msg2 = ("\n *** WARNING: *** \n You are using observer-frame "
                    "fluxes. Eiso requires rest-frame fluxes. Please use "
                    "redshift if computing fluxes for Eiso. \n")
            print(msg2)

        if isinstance(parameters, dict) is False:
            raise Exception("parameters must be a dictionary.")

        model_parts = model.split('+')
        totalFlux = []
        # Each model component's flux is computed separately and then summed. 
        for m in model_parts:
            # parlst - list of parameter values in order that funciton needs them to be. 
            pars = parameters[m].values()
            func = 'lambda x,values: x * %s(x, *values)'%(m) # x - energy
            func = eval(func)
            Flux = integrate.quad(func, emin, emax, args=(pars), limit=100)[0]
            if photonflux is True:
                totalFlux.append(Flux)
            else:
                nrgFlux = Flux*keVtoerg  # photon flux to energy flux conversion.
                totalFlux.append(nrgFlux)
        return np.sum(totalFlux)

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


        
