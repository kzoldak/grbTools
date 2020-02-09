from __future__ import division
from collections import OrderedDict
from scipy import integrate
#from math import pi
import numpy as np

#from .models import *

#from grbTools.Cosmology.luminositydistance import LumDist


def calc_flux(model, parameters, program=None, emin=None, emax=None, 
              redshift=None, photonflux=False):
    if program is not None:
        program = program.lower().replace(' ','')
        if (program == 'xspec') or (program == 'bxa'):
            # They use identical functions, however prop of uncertainty is different. 
            pass
        elif program == 'rmfit':
            pass
        else:
            # Otherwise, use general model functions.
            from general.models import *
            model_options = ["lpow", "copl", "band", "sbpl", "bbody"]
            if model not in model_options:
                raise Exception("'model' must be one of following: %r"%model_options)
    
        

    modelOptions = 'BAND SBPL COPL BBODY LPOW'.split()

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
        raise Exception("parameters must be a dictionary.")

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
        if photonflux is True:
            totalFlux.append(Flux)
        else:
            nrgFlux = Flux*keVtoerg
            totalFlux.append(nrgFlux)
    return np.sum(totalFlux)


#def Calc_Flux(model, parameters, program=None, emin=None, emax=None, 
#              redshift=None, photonflux=False):
#    """
#    Calc_Flux(model, parameters, emin, emax, redshift=None, photonflux=False)
#
#    Parameters:
#    -----------
#    model : str, model name. Full name with + between additive models.
#    parameters : ordered dict, see below.
#    program : str. Model program to grab certain functional forms. 
#              options: 'rmfit', 'xspec', 'bxa', or None.
#              ** Default is None, which uses a generalized set of functions. 
#              You almost NEVER need to specify program. The generalized form 
#              is better to use as it is more flexible. It allows you to 
#              specify which characteristic energy you have, etc. 
#              The only time you should specify a program is when you are 
#              computing fluxes for the first time (as read for a FITS file), 
#              or you propagating uncertainty (which again, you need the 
#              data from a FITS file).
#    emin : float, min energy integrand.
#    emax : float, max energy integrand.
#    redshift : float, redshift. Default is None. 
#                If redshift=None, then an un k-corrected energy flux 
#                or photon flux will be returned. 
#    photonflux : True or False. 
#                If photonflux=False, energy flux will be returned.
#                If photonflux=True, a photon flux will be returned.
#                ** Default is photonflux=False.
#
#    
#    This function calculates energy flux from an xspec model fit. 
#    We build a dictionary of the model and parameters and then pass it 
#    to Calc_Flux. This particular Calc_Flux function was set up for XSPEC 
#    models only. 
#
#    pars dictionary
#    ---------------
#    pars is an ordered dictionary, meaning the keys and values stay in the 
#      order they were entered. This is critical for this program because the 
#      Calc_Flux function was set up so that a list of parameters is passed 
#      to the model functions, not a dictionary. 
#
#      For example:
#      grbm(energy, alpha, beta, tem, norm) 
#      is the model set up for the band function. 
#
#      lpow(plIndex, norm) 
#      is the model set up for the power-law function. 
#
#      Within the pars dictionary, the parameters must be set up as so:
#          pars = OrderedDict()
#
#          pars['grbm'] = OrderedDict()
#          pars['grbm']['alpha__1']= -1.0
#          pars['grbm']['beta__2']= -2.3
#          pars['grbm']['tem__3']= 300.0
#          pars['grbm']['norm__4']= 0.01
#
#          pars['lpow'] = OrderedDict()
#          pars['lpow']['plIndex__5']= -1.0
#          pars['lpow']['norm__6']= 0.001
#
#      Notice the order of the parameters. 
#      If we do 
#          pars['grbm']['beta__2']= -2.3
#          pars['grbm']['alpha__1']= -1.0
#          pars['grbm']['tem__3']= 300.0
#          pars['grbm']['norm__4']= 0.01
#      instead, then the band function reads the -2.3 as alpha, not beta. 
#      Thus, having the parameter names defined within the dictionary is 
#      not there to help the model functions (being integrated) know the 
#      proper parameter values, but to keep things organized on our end. 
#      
#      You could do the following, and it wouldn't matter:
#          pars['grbm']['alp']= -2.3
#          pars['grbm']['bet']= -1.0
#          pars['grbm']['ech']= 300.0
#          pars['grbm']['amp']= 0.01
#      But we don't recommend doing this. 
#
#
#      Setting up models
#      ------------------
#      Regardless of whether or not you have additive models, you should 
#      always define the model name as the first level of keys in the nested 
#      dictionary. 
#      For example, if we only have the band function results, we would still 
#      do:
#            pars['grbm'] = OrderedDict()
#            pars['grbm']['alpha__1'] = -1.0
#            pars['grbm']['beta__2'] = -2.3
#            pars['grbm']['tem__3'] = 300.0
#            pars['grbm']['norm__4'] = 0.01
#
#
#    """
#    modelOptions = 'BAND SBPL COPL BBODY LPOW'.split()
#
#    keVtoerg = 1.60217657E-9
#
#    if redshift is not None:
#        emin = emin/(1.+redshift)
#        emax = emax/(1.+redshift)
#    else:
#        emin = emin
#        emax = emax
#        msg2 = ("\n *** WARNING: *** \n You are using observer-frame "
#                "fluxes. Eiso requires rest-frame fluxes. Please use "
#                "redshift if computing fluxes for Eiso. \n")
#        print(msg2)
#
#    if isinstance(parameters, dict) is False:
#        raise Exception("parameters must be a dictionary.")
#
#    #from grbTools.Models.XSPEC.models import *
#    # from grbTools.Models.XSPEC.models import grbm,sbpl,cutoffpl,bbody,lpow,powerlaw
#    # Functions are specific to XSPEC and take only those parameters. 
#    # grbm(energy, alpha, beta, tem, norm)
#    # sbpl(energy, alpha, beta, ebreak, norm)   # my function
#    # cutoffpl(energy, PhoIndex, HighECut, norm)
#    # bbody(energy, kT, norm)
#    # lpow(energy, plIndex, norm)               # my function
#    # powerlaw(energy, PhoIndex, norm)
#    modelParts = model.split('+')
#    totalFlux = []
#    for modpart in modelParts:
#        # x - energy
#        # parlst - list of parameter values in order that funciton 
#        #   needs them to be. 
#        parlst = parameters[modpart].values()
#        func = eval('lambda x,values: x * %s(x, *values)'%(modpart))
#        Flux = integrate.quad(func, emin, emax, args=(parlst), limit=100)[0]
#        if photonflux is True:
#            totalFlux.append(Flux)
#        else:
#            nrgFlux = Flux*keVtoerg
#            totalFlux.append(nrgFlux)
#    return np.sum(totalFlux)




def Calc_Fluence(model, parameters, emin, emax, duration, redshift, program):
    '''
    Calc_Fluence(model, parameters, emin, emax, duration, redshift, program)

    Parameters:
    -----------
    model: str, model name.
    parameters: ordered dict, parameter values. 
    emin: float, min energy of integration.
    emax: float, max energy of integration. 
    duration: float, duration of burst.
    redshift: float, redshift of burst.
    program: str, 'xspec', 'bxa', 'rmfit', or 'general'. 

    '''
    if (program == 'xspec') or (program == 'bxa'):
        # Somehow this is allowed, but we couldn't do this 
        #   for writing Calc_Flux and importing the model functions. 
        from grbTools.XSPEC.flux import Calc_Flux
    flux = Calc_Flux(model=model, 
                     parameters=parameters, 
                     emin=emin, 
                     emax=emax, 
                     redshift=redshift, 
                     photonflux=False)
    fluence     = flux * duration
    return fluence 


def Calc_Eiso(model, parameters, emin, emax, duration, redshift, program):
    '''
    Calc_Eiso(model, parameters, emin, emax, duration, redshift, program)

    Parameters:
    -----------
    model: str, model name.
    parameters: ordered dict, parameter values. 
    emin: float, min energy of integration.
    emax: float, max energy of integration. 
    duration: float, duration of burst.
    redshift: float, redshift of burst.
    program: str, 'xspec', 'bxa', 'rmfit', or 'general'. 

    Notes:
    ------
    Luminosity Distance Calculation can be found here:
    from grbTools.Cosmology.luminositydistance import LumDist
    Import LumDist and read its documentation to learn more on the 
    cosmology we use. 

    '''
    if (program == 'xspec') or (program == 'bxa'):
        # Somehow this is allowed, but we couldn't do this 
        #   for writing Calc_Flux and importing the model functions. 
        from grbTools.XSPEC.flux import Calc_Flux
    flux = Calc_Flux(model=model, 
                     parameters=parameters, 
                     emin=emin, 
                     emax=emax, 
                     redshift=redshift, 
                     photonflux=False)
    fluence     = flux * duration
    dL          = LumDist(redshift)
    eiso        = ((4.0 * pi * (dL**2))/(1.+redshift)) * fluence
    return eiso

