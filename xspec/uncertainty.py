"""
The idea behind propagating uncertainty:
- get string equations of the partials wrt each param. 
- take each one, evaluate the string and save it within a function.
i.e., 
pd_wrt_alpha_Func = lambda energy: eval(pd_wrt_alpha)
where pd_wrt_alpha is the partial derivative of 'grbm' wrt alpha, 
in string format. Not yet evaluated with real numbers, in numerical 
form still. 
- Then integrate that function across the desired energy range.
i.e.,  
Partials[par] = integrate.quad(pd_wrt_alpha_Func, emin, emax)[0]
- We save each one of these in a Partials dictionary. Partials should 
have a key for each model parameter. 
- Then we make a partial derivatives matrix:
  PDMAT = numpy.matrix(Partials.values())
- Take the fit covariance matrix:
COVARMAT = numpy.matrix(COVARMAT)
- Uncertainty in the model flux is then:

flux_variance = PDMAT * COVARMAT * PDMAT.T()
or 
flux_variance = PDMAT * COVARMAT * PDMAT.transpose()
then:
flux_uncertainty = float(numpy.sqrt(flux_variance)) * keVtoerg

For the Band funtion ('grbm'), there is a lower energy part and an 
upper energy part, which is due to the conditional. This complicates 
its propagation of error. You have to integrate both parts separately. 
The lower energy part is integrated from emin to engcondition and the 
upper energy part is integrated from engcondition to emax. 
engcondition stands for that energy condention that separates the 
lower and higher energy parts. 
That econdition is   ((alpha__1 - beta__2) * tem__3), aka epeak. 
The so called epeak energy is band's conditional. 
"""

from __future__ import division
from collections import OrderedDict
from scipy import integrate
import numpy as np
from math import pi, log10, log, exp, sqrt, atanh

from grbTools.Cosmology.luminositydistance import LumDist

from .modelpartialderivatives import *


def delete_pars(parameters):
    """
    Parameters:
    ------------
    parameters: dict of parameters and values. 
                parameters.keys() are the parameter names;
                'alpha__1', 'beta__2', etc.

    Clears all global parameter assignments.
    del alpha__1
    del beta__2
    del tem__3 
    etc ...

    """
    for key in parameters.keys():
        del globals()[key]


def Calc_Uncertainty(model, parameters, emin, emax, redshift, 
                     duration, covarmat, energetic):
    """
    Parameters:
    -----------
    model: str, model name. Options are:
            'grbm', 'cutoffpl', 'sbpl', 'lpow', 'bbody', 
            'grbm+lpow', 'cutoffpl+lpow', 'sbpl+lpow', 'bbody+lpow',
            'grbm+bbody', 'cutoffpl+bbody', 'sbpl+bbody',
            'grbm+bbody+lpow', 'cutoffpl+bbody+lpow', 'sbpl+bbody+lpow'
    parameters: ordered dict, holds parameter names associated with 
                the respective model. See bottom for examples. 
    emin: float, minimum rest-frame energy for integration, 
            typically 1 or 10 keV. Will be corrected for redshift. 
    emax: float, maximum rest-frameenergy for integration, 
            typically 10 MeV. Will be corrected for redshift. 
    redshift: float, redshift of GRB.
    duration: float, duration of GRB. Typically T90.
    covarmat: np.ndarray or np.matrix, 
                square covariance matrix from model fit. 
                Needs to be read from a FITS file. 
    energetic: str, 'flux', 'fluence', 'eiso'. 
                Determines which uncertainty to calculate. 
    
    Notes:
    ------
    This function can calculate flux, fluence, or eiso uncertainty. 
    It requires a covariance matrix to propagate the parameter 
    uncertainties to the model function uncertainties (and thus flux).
    Covarmats can be found in the FITS files holding their 
    XSPEC fitting results. 
    
    Future Work:
    ------------
    Not set up to handle bbody with low kTs. 
    Need to add a limit on the bbody fn's part deriv wrt kT. 
    emax = (354.8*kT)*(1+z) 
    anything above this energy will cause a math domain error. 
    This is true for all models with bbody as the pd wrt kT is 
    the same in all. 

    Examples:
    ----------
    burst = 'bn100728095'  # just for a reminder

    model = 'grbm'
    det = 'L'
    version = '-01-'
    dur = 165.378
    z = 1.567
    energetic = 'flux'


    # *~*~*~ Find appropriate files. 
    detDir = ('GBMwLAT' if 'L' in det else 'GBM')
    paramfile = ('/Users/KimiZ/GRBs2/analysis/LAT/%s/PYXSPEC/'
                 '%s/%s/xspec_fitresults_%s_%s_%s_.json'%(burst, detDir, 
                                                          model, model, 
                                                          version, det))
    covarfile = ('/Users/KimiZ/GRBs2/analysis/LAT/%s/PYXSPEC/'
                 '%s/%s/xspec_fitresults_%s_%s_%s_.fit'%(burst, detDir, 
                                                         model, model, 
                                                         version, det))
    # *~*~*~ Covaraince Matrix
    f1 = pyfits.open(covarfile)
    COVARMAT = np.asarray(f1[2].data['COVARMAT'][0])
    nPars = int(np.sqrt(len(COVARMAT)))
    COVARMAT = COVARMAT.reshape((nPars, nPars))
    #COVARMAT = np.asmatrix(COVARMAT)
    f1.close()
    del f1

    # *~*~*~ parameters dict -- object_pairs_hook=OrderedDict keeps it in order. 
    parameters = json.load(open(paramfile, 'r'), 
                           object_pairs_hook=OrderedDict)[0] 
    # Keep only the first value, which is param value
    for key in parameters.keys():
        parameters[key] = parameters[key][0]

    # *~*~*~ Calculate Uncertainty
    unc = Calc_Uncertainty(model=model, 
                           parameters=parameters, 
                           emin=10., 
                           emax=10000., 
                           redshift=z, 
                           duration=dur,
                           covarmat=COVARMAT, 
                           energetic=energetic)
    print(unc)





    parameters dict for 'sbpl+lpow'
    parameters = OrderedDict([
                    ('alpha__1', -0.8802654928774287),
                    ('beta__2', -2.73689692260396),
                    ('ebreak__3', 222.52902967845034),
                    ('norm__4', 0.011086817968969171),
                    ('plIndex__5', -0.5676834449617059),
                    ('norm__6', 2.980596437191683e-10)
                    ])

    parameters dict for 'grbm'
    parameters = OrderedDict([
                            ('alpha__1', -0.6925311809256781), 
                            ('beta__2', -2.8152292944897925), 
                            ('tem__3', 227.9801827416458), 
                            ('norm__4', 0.017064672462665017)
                            ])

    parameters dict for 'grbm+bbody+lpow'
    parameters = OrderedDict([
                            ('alpha__1', -0.9114572457472409), 
                            ('beta__2', -3.2845939589437663), 
                            ('tem__3', 365.17620202492884), 
                            ('norm__4', 0.011349424720842686), 
                            ('kT__5', 41.930243292714316), 
                            ('norm__6', 0.9741054686488074), 
                            ('plIndex__7', -0.8645465645499966), 
                            ('norm__8', 7.850109270330981e-09)
                            ])

    """
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
    DL = LumDist(redshift) # Luminosity Distance for Eiso

    if isinstance(parameters, dict) is False:
        raise Exception, "parameters must be a dictionary."
    # Dictionary holding str eqns for the partial derivatives of the 
    #   model function wrt each of its parameters. 
    #   They should already be in order. 
    Partials = eval('%s()'%model.replace('+', '_'))
    # Instantiate parameter variables (i.e., alpha__1 = -1.234, beta__2 = -2.34, ...)
    for key in parameters.keys():
        globals()[key] = parameters[key]
    
    # Partials_funcs: Dictionary holding the str fns in Partials dict, but 
    #  in function form (i.e., evaluated), so they are ready to pass to 
    #  integrate.quad()
    Partials_funcs = OrderedDict()
    nrgFluxes = OrderedDict()
    nrgFluences = OrderedDict()
    Eisos = OrderedDict()
    # Save the Partial Derivatives as functions of energy to be integrated.
    for key in Partials.keys():
        if 'grbm' in model:
            # band fn has ower and upper fns 
            Partials_funcs[key] = [lambda energy: eval(Partials[key][0]), 
                                   lambda energy: eval(Partials[key][1])]
        else:
            Partials_funcs[key] = lambda energy: eval(Partials[key])
    for key in Partials_funcs.keys():
        if 'grbm' in model:
            # band fn has ower and upper fns 
            func = lambda energy: np.piecewise(energy, 
                                               [energy < (alpha__1-beta__2)*tem__3, 
                                                energy >= (alpha__1-beta__2)*tem__3], 
                                               [Partials_funcs[key][0], 
                                                Partials_funcs[key][1]])
        else:
            func = Partials_funcs[key]
        nrgFlux = integrate.quad(func, emin, emax, limit=100)[0] * keVtoerg
        nrgFluxes[key] = nrgFlux               # energy flux
        nrgFluences[key] = nrgFlux * duration  # energy fluence
        Eisos[key] = nrgFlux * ((4.*pi*(DL**2)*duration)/(1.+redshift)) # eiso
        
    if energetic.lower() == 'flux':
        pdmatrix = np.asmatrix([nrgFluxes.values()])
    elif energetic.lower() == 'fluence':
        pdmatrix = np.asmatrix([nrgFluences.values()])
    elif energetic.lower() == 'eiso':
        pdmatrix = np.asmatrix([Eisos.values()])
    # pdmatrix is the square matrix holding the evaluated partial derivatives
    #  of the fn wrt each parameter. 
    # covarmat is the square covariance matrix from the XSPEC fit. 
    #   CALCULATING UNCERTAINTY (PROPAGATION OF ERROR)
    variance = pdmatrix * covarmat * pdmatrix.T
    uncertainty = float(np.sqrt(variance))
    delete_pars(parameters) # clears global definition of alpha__1, beta__2, ...
    return uncertainty


def Calc_Flux_Uncertainty(model, parameters, emin, emax, redshift, 
                                duration, covarmat):
    """
    see Calc_Uncertainty.__doc__
    """
    return Calc_Uncertainty(model, parameters, emin, emax, redshift, 
                            duration, covarmat, energetic='flux')


def Calc_Fluence_Uncertainty(model, parameters, emin, emax, redshift, 
                                duration, covarmat):
    """
    see Calc_Uncertainty.__doc__
    """
    return Calc_Uncertainty(model, parameters, emin, emax, redshift, 
                            duration, covarmat, energetic='fluence')


def Calc_Eiso_Uncertainty(model, parameters, emin, emax, redshift, 
                                duration, covarmat):
    """
    see Calc_Uncertainty.__doc__
    """
    return Calc_Uncertainty(model, parameters, emin, emax, redshift, 
                            duration, covarmat, energetic='eiso')






def Calc_Epeak_Uncertainty(model, parameters, covarmat):
    """
    Calculates Epeak Uncertainty.

    Parameters:
    -----------
    model: str, model. If you have 'grbm' with 'tem', this gives you
                uncertainty in epeak. Need to compute epeak elsewhere. 

    parameters: ordered dict, holds parameter names associated with 
                the respective model. See bottom for examples. 
    covarmat: np.ndarray or np.matrix, 
                square covariance matrix from model fit. 
                Needs to be read from a FITS file. 

    Notes:
    ------
    The function for calculating epeak from ebreak or tem does not 
    include the model amplitude term (i.e., norm__4 parameter). 
    Therefore it doesn't make sense to include a partial derivative of 
    the fn wrt this parameter. HOWEVER, we did this, which is 0, so that 
    we can leave the covariance matrix from the fit in its present form. 
    If we removed the pd wrt norm__4, then we'd have to remove all the 
    norm__4 terms in the covariance matrix (which is all of row 3 
    and col 3 for sbpl and grbm fits - starting at 0). Without the 
    norm__4 in the partial derivative matrix, the shapes between the 
    pd matrix and covariance matrix will be different.  
    We tested it, and if we remove the norm__4 terms from both the 
    pd matrix and the covariance matrix, we will get the same thing if 
    we just leave the terms in both. 


    """

    if isinstance(parameters, dict) is False:
        raise Exception, "parameters must be a dictionary."
    # Dictionary holding str eqns for the partial derivatives of the 
    #   model function wrt each of its parameters. 
    #   They should already be in order. 
    Partials = eval("epeak_partials('%s')"%model)

    # Instantiate only appropriate parameters. 
    #  replaced parameteres.keys() with Partials.keys()
    for key in Partials.keys():
        globals()[key] = parameters[key]

    # Evaluate partials, use same dict. No need to change it. 
    for key in Partials.keys():
        Partials[key] = eval(Partials[key])  # holds real numbers now. 
    pdmatrix = np.asmatrix([Partials.values()])  # make matrix.
    # pdmatrix is the square matrix holding the evaluated partial derivatives.
    # covarmat is the square covariance matrix from the XSPEC fit. 
    #   CALCULATING UNCERTAINTY (PROPAGATION OF ERROR)
    variance = pdmatrix * covarmat * pdmatrix.T
    uncertainty = float(np.sqrt(variance))

    if model == 'grbm':
        epeak = (alpha__1 + 2.0) * tem__3
    elif model == 'cutoffpl':
        epeak = (PhoIndex__1 + 2.0) * HighECut__2
    elif model == 'sbpl':
        epeak = ebreak__3*(10**(0.3*(atanh((alpha__1+beta__2+4.0) \
                        /(alpha__1-beta__2)))))
    else:
        raise Exception, "Don't recognize model."
    delete_pars(parameters) # clears global definition of alpha__1, beta__2, ...
    return epeak, uncertainty



def Calc_Epeak_Uncertainty_Assym(model, parameters, covarmat_lo, covarmat_up):
    """
    Calculates Epeak Uncertainty.

    Parameters:
    -----------
    model: str, model. If you have 'grbm' with 'tem', this gives you
                uncertainty in epeak. Need to compute epeak elsewhere. 

    parameters: ordered dict, holds parameter names associated with 
                the respective model. See bottom for examples. 
    covarmat_lo: np.ndarray or np.matrix, 
                 square covariance matrix from model fit, but updated 
                 so that the variance terms (diag terms) are replaced 
                 with the square of the lower margin of error on the 
                 parameter. This needs to be done before passing 
                 covariance matrices to this function. 
    covarmat_up: np.ndarray or np.matrix, 
                 square covariance matrix from model fit, but updated 
                 so that the variance terms (diag terms) are replaced 
                 with the square of the upper margin of error on the 
                 parameter. This needs to be done before passing 
                 covariance matrices to this function. 

    BEFORE PASSING COVARIANCE MATRICES TO THIS FUNCTION DO THE FOLLOWING:
    COVARMAT_LO = COVARMAT.copy()
    COVARMAT_UP = COVARMAT.copy()
    for i,key in enumerate(parameters.keys()):
        COVARMAT_LO[i][i] = parameters[key][1]**2
        COVARMAT_UP[i][i] = parameters[key][2]**2

    where parameters[parname] = [value, moe_lo, moe_up]


    Notes:
    ------
    ** WE DON'T RECOMMEND USING THIS FUNCTION. NOT ENTIRLEY SURE YOU 
        CAN DO THIS. **

    The function for calculating epeak from ebreak or tem does not 
    include the model amplitude term (i.e., norm__4 parameter). 
    Therefore it doesn't make sense to include a partial derivative of 
    the fn wrt this parameter. HOWEVER, we did this, which is 0, so that 
    we can leave the covariance matrix from the fit in its present form. 
    If we removed the pd wrt norm__4, then we'd have to remove all the 
    norm__4 terms in the covariance matrix (which is all of row 3 
    and col 3 for sbpl and grbm fits - starting at 0). Without the 
    norm__4 in the partial derivative matrix, the shapes between the 
    pd matrix and covariance matrix will be different.  
    We tested it, and if we remove the norm__4 terms from both the 
    pd matrix and the covariance matrix, we will get the same thing if 
    we just leave the terms in both. 

    """

    if isinstance(parameters, dict) is False:
        raise Exception, "parameters must be a dictionary."
    # Dictionary holding str eqns for the partial derivatives of the 
    #   model function wrt each of its parameters. 
    #   They should already be in order. 
    Partials = eval("epeak_partials('%s')"%model)

    # Instantiate only appropriate parameters. 
    #  replaced parameteres.keys() with Partials.keys()
    for key in Partials.keys():
        globals()[key] = parameters[key][0]

    # Evaluate partials, use same dict. No need to change it. 
    for key in Partials.keys():
        Partials[key] = eval(Partials[key])  # holds real numbers now. 
    pdmatrix = np.asmatrix([Partials.values()])  # make matrix.
    # pdmatrix is the square matrix holding the evaluated partial derivatives.
    # covarmat is the square covariance matrix from the XSPEC fit. 
    #   CALCULATING UNCERTAINTY (PROPAGATION OF ERROR)
    variance_lo = pdmatrix * covarmat_lo * pdmatrix.T
    uncertainty_lo = float(np.sqrt(variance_lo))

    variance_up = pdmatrix * covarmat_up * pdmatrix.T
    uncertainty_up = float(np.sqrt(variance_up))

    if model == 'grbm':
        epeak = (alpha__1 + 2.0) * tem__3
    elif model == 'cutoffpl':
        epeak = (PhoIndex__1 + 2.0) * HighECut__2
    elif model == 'sbpl':
        epeak = ebreak__3*(10**(0.3*(atanh((alpha__1+beta__2+4.0)/(alpha__1-beta__2)))))
    else:
        raise Exception, "Don't recognize model."
    delete_pars(parameters) # clears global definition of alpha__1, beta__2, ...
    return epeak, uncertainty_lo, uncertainty_up



