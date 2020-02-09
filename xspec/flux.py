from __future__ import division
from collections import OrderedDict
from scipy import integrate
from math import pi
import numpy as np

from .models import *


def Calc_Flux(model, parameters, emin, emax, redshift=None, photonflux=False):
    """
    Calc_Flux(model, parameters, emin, emax, redshift=None, photonflux=False)

    Parameters:
    -----------
    model: str, model name. Full name with + between additive models.
    parameters: ordered dict, see below.
    emin: float, min energy integrand.
    emax: float, max energy integrand.
    redshift: float, redshift.
    photonflux: True or False, default is False. True returns photon fluxes 
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
    modelOptions = 'grbm sbpl cutoffpl bbody lpow powerlaw'.split()

    keVtoerg = 1.60217657E-9

    if redshift is not None:
        emin = emin/(1.+redshift)
        emax = emax/(1.+redshift)
    else:
        emin = emin
        emax = emax
        msg2 = ("\n *** WARNING: *** \n You are using observer-frame "
                "fluxes. Eiso requires rest-frame fluxes. Please use "
                "redshift if computing fluxes for Eiso. \n")
        print(msg2)

    if isinstance(parameters, dict) is False:
        raise Exception, "parameters must be a dictionary."

    #from grbTools.Models.XSPEC.models import *
    # from grbTools.Models.XSPEC.models import grbm,sbpl,cutoffpl,bbody,lpow,powerlaw
    # Functions are specific to XSPEC and take only those parameters. 
    # grbm(energy, alpha, beta, tem, norm)
    # sbpl(energy, alpha, beta, ebreak, norm)   # my function
    # cutoffpl(energy, PhoIndex, HighECut, norm)
    # bbody(energy, kT, norm)
    # lpow(energy, plIndex, norm)               # my function
    # powerlaw(energy, PhoIndex, norm)
    modelParts = model.split('+')
    totalFlux = []
    for modpart in modelParts:
        # x - energy
        # parlst - list of parameter values in order that funciton 
        #   needs them to be. 
        parlst = parameters[modpart].values()
        func = eval('lambda x,values: x * %s(x, *values)'%(modpart))
        Flux = integrate.quad(func, emin, emax, args=(parlst), limit=100)[0]
        nrgFlux = Flux*keVtoerg
        if photonflux is True:
            totalFlux.append(Flux)
        else:
            totalFlux.append(nrgFlux)
    return np.sum(totalFlux)


