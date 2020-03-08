"""
This is for GENERAL USE ONLY! 
Most of the time you want to use these!
Any time you are calculating flux, fluence, or eiso for the first time 
and reading results from a program's FITS files, you will want the functions 
within the xspec and rmfit directories!



The model functions in this file are flexible and were written for GENERIC 
use, not to be used with FITS file results!!!!

Most of the time, you will use these functions! If you need the ones 
tailored exactly for XSPEC or RMFIT, use the files within the xspec and rmfit 
folders. There are calc_flux, calc_fluence, and calc_eiso functions for those, 
and the models are specifically designed to be used with their respecitve 
FITS files. 



When to use these functions:
----------------------------
- When you have results NOT coming directly from a FITS file. 
- When you do not have a covariance matrix for uncertainty propagation.
- When you want a quick flux, fluence, or Eiso; or a quick plot of the model.


When NOT to use these functions:
--------------------------------
- When you are calculating flux, fluence, or Eiso for the first time.
- When your results are coming directly from an XSPEC or RMFIT produced 
  FITS file. 
- When you need accurate propagation of uncertainty. 



Notes on these functions:
-------------------------
Model options are:
    band, sbpl, copl, bbody, and lpow 
Any time the user enters 'grbm' instead of 'band', the band function 
is called and the user is warned that these are general functions and that 
'band' is being used. Same thing for copl vs cutoffpl. 



When you have a set of parameters and you want to calculate a flux, 
fluence, or Eiso, or make a plot of the model in nuFnu space, then 

They were written with flexibility. 
Function normalizations chan be changed, so can characteristic energy. 

These functions were designed with XSPEC in mind, but can be adapted to match 
RMFIT's functional form. 

The SBPL model does not exist within XPSEC, so I wrote it identical to RMFIT's.

The three main models: BAND, SBPL, and COPL  were designed to take different 
characteristic energies, you just have to specify which one you are passing.
   For example, you can pass BAND 'epeak', 'tem', 'E0', or 'ebreak'. 
   'tem' and 'E0' mean the same thing. Default is 'E0'. 
   COPL can take 'epeak', 'ebreak', 'E0', or 'HighECut'.

"""
from __future__ import division # just incase this is running on Python 2.7
from collections import OrderedDict

import warnings

from scipy import integrate
from math import pi
import numpy as np

from numpy import exp, log, log10
from math import atanh

from .models import *

from grbTools.errors.errorhandling import RedshiftWarning, NotADictError

from grbTools.cosmology.luminositydistance import LumDist  # basic version



__all__ = ["calc_flux", "calc_fluence", "calc_eiso", "parameter_dict_example"]



# class RedshiftWarning(UserWarning):
#     # Warn user that redshift is not being used.
#     pass

# class NotADictError(Exception):
#     # Must use a dict.
#     pass




# *********************   FLUX FUNCTION   *****************

def calc_flux(modelname, parameters, emin, emax, redshift=None, 
    photonflux=False):
    """
    This is for GENERAL USE ONLY! 
    Most of the time you want to use these!
    Any time you are calculating flux, fluence, or eiso for the first time 
    and reading results from a program's FITS files, you will want the functions 
    within the xspec and rmfit directories!


    Parameters:
    -----------
    modelname : str, model name. Full name with + between additive models.
    parameters : ordered dict, see below.
    emin : float, min energy integrand.
    emax : float, max energy integrand.
    redshift : float, redshift. Default is None. 
                If redshift=None, then an un k-corrected energy flux 
                or photon flux will be returned. 
    photonflux : True or False. 
                If photonflux=False, energy flux will be returned.
                If photonflux=True, a photon flux will be returned.
                ** Default is photonflux=False.


    parameters dict.
    ----------------
    For each model component in the modelname, the fluxes are calculated 
    separately and them summed together. 
    For example, if you modelname is 'grbm+bbody+lpow', then the flux for 
    'grbm' is calcualted first and stroed in a Flux list, 
    then 'bbody' flux is calculated second and stored in the Flux list,
    then 'lpow' flux is calculated and stored in the Flux list. 
    At the end, the Total Flux is found by summing the values in the list 
    together. 
    TotalFlux = sum(Fluxlist)

    Hence the reason for passing the parameters for each model component 
    in a sub-dict within the main parameters dict. 
       E.g., 
    parameters = {'grbm': {'alpha': -1.23, 'beta':-2.56, ...}, 
                  'bbody': {'kT':54.3, ...}, 
                  'lpow': {'plIndex':-1.49, ...}
                  }
    we will provide more accurate examples of these in a moment. 

    It is recommended that you use OrderedDict from collections

    from collections import OrderedDict

    parameters = OrderedDict()

    parameters['grbm'] = OrderedDict()   # sub-dict for 'grbm'

    parameters['grbm']['alpha'] = -1.0
    parameters['grbm']['beta'] = -2.5
    parameters['grbm']['enterm'] = 500.0
    parameters['grbm']['norm'] = 0.01
    parameters['grbm']['entype'] = 'epeak'  # This is how you determine which char. eng term.

    parameters['bbody'] = OrderedDict() 
    parameters['bbody']['kT'] = 54.1
    parameters['bbody']['norm'] = 0.1

    parameters['lpow'] = OrderedDict()   # sub-dict for lpow
    parameters['lpow']['plIndex']= -1.0
    parameters['lpow']['norm']= 0.001

    The exact parameter names are insignificant, just as long as they 
    are in the correct order. They are flattened into a list and then 
    passed to the functions. 

    parameters['grbm']['entype']


    !!!!!!!! WARNING !!!!!!!!
    The 

    """
    if isinstance(parameters, dict) is False:
        raise NotADictError("'parameters' must be a dictionary.")

    for m in modelname.split('+'):
        if m not in parameters.keys():
            msg = (" Use nested dicts for passing model parameters via "
                   "the 'parameters' attribute. "
                   " Run \n\tparameter_dict_example('%s')"
                   " for an example."%modelname)
            raise NotADictError(msg)

        keVtoerg = 1.60217657E-9
        if redshift is not None:
            emin = emin/(1.+redshift)
            emax = emax/(1.+redshift)
        else:
            msg2 = ("\n *** WARNING: *** \n You are using observer-frame "
                    "fluxes. Eiso requires rest-frame fluxes. Please use "
                    "redshift if computing fluxes for Eiso. \n")
            print(msg2)

        # Each model component's flux is computed separately and then summed. 
        model_parts = modelname.split('+')
        totalFlux = []
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


    modelOptions = ["lpow", "copl", "band", "grbm", "sbpl", "bbody"]
    keVtoerg = 1.60217657E-9

    if redshift is not None:
        emin = emin/(1.+redshift)
        emax = emax/(1.+redshift)
    else:
        #warnings.warn(msg, UserWarning, stacklevel=1)
        # msg = ("\n *** WARNING: *** \n You are using observer-frame "
        #         "fluxes. Eiso requires rest-frame fluxes. Please use "
        #         "redshift if computing fluxes for Eiso. \n")
        msg = ("\nNo redshift provided, calculating observer-frame flux. "
               "If calculating for Eiso energy, pass a redshift. Eiso "
               "is a rest-frame energetic.")
        warnings.warn(msg, RedshiftWarning, stacklevel=1)
        #print(msg)


    # totalFlux = []
    # for mod in modelname.split('+'):  # flux calcualted for each component separately.
    #     # parlst - list of parameter values in order that funciton 
    #     #   needs them to be. 
    #     if mod not in parameters.keys():
    #         msg = "Make sure you use nested dict for passing 'parameters'. "
    #         raise Exception(msg)

    #     if not isinstance(parameters, OrderedDict):
    #         msg = ("Please use an OrderedDict. If not, make sure you enter parameters in the correct order. ")
    #         warnings.warn(msg, UserWarning, stacklevel=1)

    #     #pars = parameters[mod].values()

    #     pars = parameters[mod]['alpha']
    #     pars = parameters[m].values()


    #     func = eval('lambda x,values: x * %s(x, *values)'%(mod))  # x - energy
    #     flux = integrate.quad(func, emin, emax, args=(pars), limit=100)[0]
    #     if photonflux is True:
    #         totalFlux.append(flux)
    #     else:
    #         nrgFlux = flux*keVtoerg
    #         totalFlux.append(nrgFlux)
    # return np.sum(totalFlux)

