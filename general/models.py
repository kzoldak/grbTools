"""
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
from __future__ import division  # just incase this is running on Python 2.7

import warnings

import numpy as np
from numpy import exp, log, log10
from math import atanh



__all__ = ["lpow", "powerlaw", "copl", "band", "grbm", "sbpl", "bbody"]



def band(energy, alpha, beta, enterm, norm, entype='epeak', enorm=100.0):
    '''

    Band function. 

    Parameters:
    -----------
    energy : float, 
                a *single* energy from an array of energies over which this 
                model's flux is integrated. 
    alpha :  float, value of alpha parameter (i.e., low-energy index). 
    beta :   float, value of beta parameter (i.e., high-energy index).    
    enterm : float, value of characteristic energy (see parameter `entype`).  
    norm :   float, value of model amplitude (often called model normalization, 
                which can be confused with parameter `enorm`).  
    entype : str, defines which characteristic energy you are passing 
                to this function in parameter 'enterm'.  
                Options are: 
                    'epeak' (default), E0', 'tem', or 'ebreak'.
                    'E0' and 'tem' are the same.  
    enorm : int or float, optional since default value is provided. 
                Energy at which the model is normalized at. 
                Default is 100 keV, but XSPEC often normalizes its functions to 
                1 keV. This one uses enorm = 100.0 keV. 

    Notes:
    ------
    This function is set up in a general form to take any characteristic energy 
    you may have. It can also handle different enorm energies. 


    Calculations:
    -------------
    tem and E0 are the same thing. Conversion between 
    E0, epeak, and ebreak:

    E0 = epeak / (2.+alpha)
    E0 = ebreak / (alpha - beta)

    epeak = E0 * (2.+alpha)
    ebreak = E0 * (alpha - beta) = (epeak/(2.+alpha)) * (alpha-beta)
    epeak = (ebreak/(alpha-beta)) * (2.+alpha)
    
    '''
    entype_options = ['E0','tem','epeak','ebreak']
    if entype not in entype_options:
        msg = 'entype must be one of following: %r'%entype_options
        raise Exception(msg)
        
    eng = energy
    Amp = norm
    a = alpha
    b = beta
    enorm = enorm

    if (entype == 'E0') or (entype == 'tem'):
        epk = enterm * (2.+a)
    elif entype == 'ebreak':
        epk = (enterm/(a-b)) * (2.+a)
    else:
        epk = enterm  # default is entype == 'epeak'

    condLO = eng < ((a-b) * epk)/(2.0+a)
    condUP = eng >= ((a-b) * epk)/(2.0+a)

    bandLO = lambda x: Amp * pow((x/enorm), a) * exp(-(2.+a)*x/epk)
    bandUP = lambda x: Amp * pow((x/enorm), b) * exp(b-a) * \
                        pow(((a-b)*epk)/(enorm*(2.+a)), a-b)
    result = np.piecewise(eng, [condLO, condUP], [bandLO, bandUP])
    return result


def grbm(energy, alpha, beta, enterm, norm, entype='epeak', enorm=100.0):
    """
    A copy of the band function.

    """
    return band(energy, alpha, beta, enterm, norm, entype=entype, enorm=enorm)


def sbpl(energy, alpha, beta, enterm, norm, entype='ebreak', 
            enorm=100.0, brkscale=0.3):
    '''

    Sbpl function. 

    Parameters:
    -----------
    energy : float, 
                a *single* energy from an array of energies over which this 
                model's flux is integrated. 
    alpha :  float, value of alpha parameter (i.e., low-energy index). 
    beta :   float, value of beta parameter (i.e., high-energy index).    
    enterm : float, value of characteristic energy (see parameter `entype`).  
    norm :   float, value of model amplitude (often called model normalization, 
                which can be confused with parameter `enorm`).  
    entype : str, defines which characteristic energy you are passing 
                to this function in parameter 'enterm'.  
                Options are: 
                    'ebreak' (default), and 'epeak'. 
    enorm : int or float, optional since default value is provided. 
                Energy at which the model is normalized at. 
                Default is 100 keV, but XSPEC often normalizes its functions to 
                1 keV. This one uses enorm = 100.0 keV. 
    brkscale : int or float, optional since default value is provided.
                Sbpl's break scale energy. The default value is 0.3, 
                which was found to be an optimal value by Kaneko et al (??). 
                They did simulations on CGRO/BATSE bursts and fit the data 
                with this model in search of a value that was unique enough 
                from Band to merit its use. 

    Notes:
    ------
    This function is set up in a general form to take any characteristic energy 
    you may have. It can also handle different enorm energies. 

    ****
    This function does not exist within XSPEC. We created a user-defined 
    Python function for Sbpl to be adopted by PyXspec for spectral fitting. 
    That function matches this one exactly, excpet that this formalism is a 
    bit more flexible. In the PyXspec version, we fixed brkscale = 0.3 and 
    we do not provide the option of another enorm. 

    We made sure this funtion also matched exactly to the RMFIT version. 


    Calculations:
    -------------
    Conversion between epeak and ebreak:
    a = alpha
    b = beta
    d = 0.3   # break scale.
    epeak   = (ebreak) * (10**(d * atanh((a+b+4)/(a-b))))
    ebreak  = (epeak) / (10**(d * atanh((a+b+4)/(a-b))))

    pow(10, d) is same as 10**d

    '''
    entype_options = ['epeak','ebreak']
    if entype not in entype_options:
        msg = 'entype must be one of following: %r'%entype_options
        raise Exception(msg)

    eng  = energy
    a    = alpha
    b    = beta
    d    = brkscale
    Amp  = norm
    enorm = enorm
    
    
    if entype == 'epeak':
        k = (enterm) / (10**(d * atanh((a+b+4.)/(a-b))))
    else:
        k = enterm  # default is to use ebreak

    p1   = (b-a)/2.
    p2   = (a+b)/2.
    p3   = log10(enorm/k)/d
    p4   = log10(eng/k)/d

    # pow(10, d) is same as 10**d
    part1 = p1 * d * log((exp(p4) + exp(-p4))/2.)
    part2 = p1 * d * log((exp(p3) + exp(-p3))/2.)
    
    result = Amp * pow(eng/enorm, p2) * pow(10, part1-part2)
    return result
# result = Amp * ((eng/enorm)**p2) * (10**((p1 * d * log((exp(p4) + exp(-p4))/2.)) \
#     - (p1 * d * log((exp(p3) + exp(-p3))/2.))))



def copl(energy, alpha, enterm, norm, entype='epeak', enorm=100.0):
    '''
    Cutoff Power-law function. This function reflect's RMFIT's version; 
        see Notes!!

    Parameters:
    -----------
    energy : float, 
                a *single* energy from an array of energies over which this 
                model's flux is integrated. 
    alpha :  float, value of alpha parameter (i.e., low-energy index). 
    enterm : float, value of characteristic energy (see parameter `entype`).  
    norm :   float, value of model amplitude (often called model normalization, 
                which can be confused with parameter `enorm`).  
    enterm : float, value of characteristic energy (see parameter `entype`).  
    norm :   float, value of model amplitude (often called model normalization, 
                which can be confused with parameter `enorm`).  
    entype : str, defines which characteristic energy you are passing 
                to this function in parameter 'enterm'.  
                Options are: 
                    'epeak' (default), E0', 'tem', or 'ebreak'.
                    'E0' and 'tem' are the same.  
    enorm : int or float, optional since default value is provided. 
                Energy at which the model is normalized at. 
                Default is 100 keV. 

                but XSPEC often normalizes its functions to 
                1 keV. This one uses enorm = 100.0 keV. 
                RMFIT:  100 keV
                XSPEC:  1 keV

    Notes:
    ------
    This function is set up in a general form to take any characteristic energy 
    you may have. It can also handle different enorm energies. 

    *** XSPEC versus RMFIT ***
    This function is setup like RMFIT's. 
    This function is normalized (enorm) to 100 keV, whereas XSPEC's version 
      is normalized to 1 keV. 
    In addition, the sign on alpha here is oposite of XSPEC's 
      plIndex; alpha and plIndex represent the same thing. 

    If you have XSPEC derived cutoffpl parameters, the conversion can be 
      done in one of the two following ways: 

        alpha = -alpha
        norm = 100^(-alpha)  # This alpha is AFTER you swap alpha's sign. 
        
        # or to do the swap simultaneously in Python,

        alpha,norm = -alpha, 100**alpha  
        # The tuple on the right side is evaluated before the new assignments. 

    The copl model is the band function in the limit that beta goes to neative 
      infinity. Which means that you can set this function up as a piecewise 
      based on two conditions (as band's function is), and make the output 
      of the second conditional (higher energies) to be zero. We didn't 
      do that here.



    Calculations:
    -------------
    tem and E0 are the same thing. Conversion between 
    E0, epeak, and ebreak:

    E0 = epeak / (2.+alpha)
    E0 = ebreak / (alpha - beta)

    epeak = E0 * (2.+alpha)
    ebreak = E0 * (alpha - beta) = (epeak/(2.+alpha)) * (alpha-beta)
    epeak = (ebreak/(alpha-beta)) * (2.+alpha)
                
    
    '''
    entype_options = ['E0','tem','epeak','ebreak']
    if entype not in entype_options:
        msg = 'entype must be one of following: %r'%entype_options
        raise Exception(msg)
        
    eng = energy
    Amp = norm
    a = alpha
    enorm = enorm

    if (entype == 'E0') or (entype == 'tem'):
        epk = enterm * (2.+a)
    elif entype == 'ebreak':
        epk = (enterm/(a-b)) * (2.+a)
    else:
        epk = enterm  # default is entype == 'epeak'

    # E0 = epeak / (2.+alpha)
    result = Amp * pow(eng/enorm, a) * exp(-(2.+a)*eng/epk)
    return result



def bbody(energy, kT, norm, program='xspec'):
    '''
    Blackbody function. This function reflect's XSPEC's version; 
        see Notes!!


    Parameters:
    -----------
    energy : float, 
                a *single* energy from an array of energies over which this 
                model's flux is integrated. 
    kT:     float, value of blackbody temperature in keV. 
    norm :   float, value of model amplitude (often called model normalization, 
                which can be confused with parameter `enorm`). 


    Notes:
    ------
    This function is not setup with an enorm parameter. 
    We fashioned it after the XSPEC version, not the RMFIT version.


    Calculations:
    -------------
    XSPEC version
    Amp * (pow(eng, 2) * 8.0525)/(pow(temp, 4) * (exp(eng/temp)-1.))

    RMFIT version
    Amp * pow(eng, 2) * pow(exp(eng/temp)-1, -1)

    The difference between the two is a factor of 
    8.0525/kT**4

    To converte between model amplitues:
    A_r = A_x * (8.0525/kT**4)
    where A_r is the RMFIT amplitude 
    and 
    A_x is XSPEC's amplitude.
    XSPEC's is normalized to 1 keV whereas RMFITs is normalized to 100 keV.

    When passing an XSPEC norm value, it should be > 0.5 and RMFIT's should be 
    on the order of < 1E-4.

    '''
    eng     = energy
    temp    = kT
    Amp     = norm
    if program == 'xspec':
        result = Amp * (pow(eng, 2) * 8.0525)/(pow(temp, 4) * (exp(eng/temp)-1.))
    elif program == 'rmfit':
        result = Amp * (pow(eng, 2)/(exp(eng/temp)-1.))
    else:
        raise Exception("'program' must be 'xspec' or 'rmfit'. ")
    return result

    #if eng <= (709.666 * kT): # to avoid exp overflow error
    #    return Amp * (((eng**2)*(8.0525)) / ((kT**4) * (exp(eng/kT)-1)))
    #else:
    #    return 0
    

def lpow(energy, alpha, norm, enorm=100.0):
    '''
    Our power-law function. 
    Same as RMFIT's. 
    Slightly different than XSPEC's: 
        - different sign on alpha 
        - normalized to 100 keV rather than 1 keV, like XSPEC's.
    
    PARAMETERS:
    ----------
    energy : float, energy element from an array of energies over which 
            the function is integrated over. 
    alpha :  float, power-law index parameter value.  
    norm :   float, normalization (aka amplitude) parameter value.  

    '''
    eng = energy
    a = alpha
    Amp = norm
    result = Amp * pow(eng/enorm, a)
    return result


def powerlaw(energy, alpha, norm, enorm=1.0):
    '''
    XSPEC power-law function. 

    Is normalized to 1 keV whereas RMFIT is to 100 keV.
    alpha sign is backwards from RMFIT's.
    
    PARAMETERS:
    ----------
    energy : float, energy element from an array of energies over which 
            the function is integrated over. 
    alpha :  float, power-law index parameter value.  
    norm :   float, normalization (aka amplitude) parameter value.  

    '''
    eng = energy
    a = alpha
    Amp = norm
    result = Amp * pow(eng/enorm, -a)
    return result





### *********************************************************************
##  BELOW:   TRYING TO MAKE FUNCTIONS ACCEPT ATTRIBUTES AND DICT FOR THE 
##              PARAMETERS. NEED TO IMPLEMENT WARNINGS EVERYWHERE. 
##          !! FUTURE WORK !!
### *********************************************************************


# def band(energy=None, alpha=None, beta=None, enterm=None, norm=None, 
#     entype=None, enorm=None, **parsDict):


#     # Read from dict.
#     if parsDict:
        
#         a = parsDict['alpha']
#         b = parsDict['beta']
#         Amp = parsDict['norm']
#         if 'enorm' in parsDict.keys():
#             enorm = parsDict['enorm']
#         else:
#             warnings.warn("Using enorm=100.0. Specify 'enorm' if this action "
#             "is not desired.", 
#             UserWarning, stacklevel=1)
#             enorm = 100.0
#     else:
#         # Assume attributed names used. 
#         if entype is None:
#             warnings.warn("Using entype='epeak'. Specify 'entype' if this "
#                 "action is not desired.", 
#                 UserWarning, stacklevel=1)
#             entype = 'epeak'

#     if enorm is None:
#         warnings.warn("Using enorm=100.0. Specify 'enorm' if this action "
#             "is not desired.", 
#             UserWarning, stacklevel=1)
#         enorm = 100.0



#     entype_options = ['E0','tem','epeak','ebreak']
#     if entype not in entype_options:
#         msg = 'entype must be one of following: %r'%entype_options
#         raise Exception(msg)
        
#     eng = energy
#     Amp = norm
#     a = alpha
#     b = beta
#     enorm = enorm

#     if (entype == 'E0') or (entype == 'tem'):
#         epk = enterm * (2.+a)
#     elif entype == 'ebreak':
#         epk = (enterm/(a-b)) * (2.+a)
#     else:
#         epk = enterm  # default is entype == 'epeak'

#     condLO = eng < ((a-b) * epk)/(2.0+a)
#     condUP = eng >= ((a-b) * epk)/(2.0+a)

#     bandLO = lambda x: Amp * pow((x/enorm), a) * exp(-(2.+a)*x/epk)
#     bandUP = lambda x: Amp * pow((x/enorm), b) * exp(b-a) * \
#                         pow(((a-b)*epk)/(enorm*(2.+a)), a-b)
#     result = np.piecewise(eng, [condLO, condUP], [bandLO, bandUP])
#     return result


# def grbm(energy, alpha, beta, enterm, norm, entype='epeak', enorm=100.0):
#     """
#     A copy of the band function.

#     """
#     return band(energy, alpha, beta, enterm, norm, entype=entype, enorm=enorm)
