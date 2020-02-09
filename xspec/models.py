'''
All of the functions in this file were designed based on the XSPEC modeling functions. 
The SBPL model does not exist within XPSEC, so we write it identical to RMFIT's.

The three main models: BAND, SBPL, and COPL  were designed to take different characteristic
energies. For example, you can pass BAND 'epeak', 'E0', or 'ebreak'. Default is 'E0'. 


'''
from __future__ import division
import numpy as np
from numpy import exp, power, log, log10, exp
from math import atanh


def grbm(energy, alpha, beta, tem, norm):
    '''

    Parameters:
    -----------
    energy : float, energy to integrate over. 
    alpha :  float, alpha (low-energy index) parameter value.  
    beta :   float, beta (high-energy index) parameter value.  
    tem : float, characteristic energy in keV. Aka E0.
            Same as the cutoffpl's HighECut.
            RMFIT uses 'epeak'. 
    norm : float, model normalization or amplitude. 

    Notes:
    -------
    This function is normalized to 100 keV and can not be changed. 

    XSPEC's 'grbm' function (i.e., the Band function).
    https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node179.html

    XSPEC's band function is the same as RMFIT's, except for:
            tem = epeak/(2.+alpha)
    In relation to ebreak: 
        tem = ebreak/(alpha-beta)

    Equation:
    ---------
    if energy < (alpha-beta)*tem:
        return norm * (((energy/100.0)**alpha) * (exp(-energy/tem)))
    else:
        return norm * (((((alpha-beta)*tem)/100.0)**(alpha-beta)) * \
            (exp(beta-alpha))*((energy/100.0)**beta))


    We set the function up as a piecewise. 

    '''
    eng     = energy
    Amp     = norm
    a       = alpha
    b       = beta
    tem     = tem

    condLo = eng < (a-b) * tem
    condUp = eng >= (a-b) * tem
    # can use power(x/100., a) instead of ((x/100.)**a)
    # Upper and Lower parts of the conditional. 
    lowerBand = lambda x: Amp * (((x/100.0)**a)*exp(-x/tem))
    upperBand = lambda x: Amp * (((((a-b)*tem)/100.0)**(a-b))*(exp(b-a))*((x/100.0)**b))
    # Return lowerBand to integrate if condLo met. 
    # Return upperBand to integrate if condUp met. 
    return np.piecewise(eng, [condLo, condUp], [lowerBand, upperBand])


def sbpl(energy, alpha, beta, ebreak, norm):
    '''
    sbpl(energy, alpha, beta, ebreak, norm)

    Parameters:
    -----------
    energy : float, energy to integrate over. 
    alpha :  float, alpha (low-energy index) parameter value.  
    beta :   float, beta (high-energy index) parameter value.    
    ebreak : float, characteristic energy in keV.   
    norm : float, model amplitude. This is normalized to 100 keV 
            in both RMFIT and XSPEC.
 

    Notes:
    ------
    This function is normalized to 100 keV and can not be changed. 

    This function does not exist within XSPEC. We created it in python. 

    We created a user-designed model model for fitting (i.e., we wrote
    this function). It is the exact same as the RMFIT's.
    We do not allow the energy break scale to vary during fitting; we 
    keep it fixed at 0.3. This function keeps it fixed at 0.3 as well. 

    If the user ever wants to change this, the parameter 'd' within 
    the code below is what needs to be changed to be a function argument. 
    For example, SBPL(energy, alpha, beta, enterm, brkscl, norm, entype='ebreak')
    where d (in func code) is the brkscl. 
    Conversion between epeak and ebreak:
    epeak   = (ebreak) * (10**(0.3 * atanh((a+b+4)/(a-b))))
    ebreak  = (epeak) / (10**(0.3 * atanh((a+b+4)/(a-b))))
          where a and b are alpha and beta

    ebreak = (enterm) / (10**(d * atanh((a+b+4.)/(a-b))))
    epeak = (ebreak) * (10**(0.3 * atanh((a+b+4)/(a-b))))
    ebreak = (epeak) / (10**(0.3 * atanh((a+b+4)/(a-b))))
    '''
    eng  = energy
    a    = alpha
    b    = beta
    k    = ebreak
    d    = 0.3  # brkscl. 
    Amp  = norm

    # To shorten the function:
    p1   = (b-a)/2.
    p2   = (a+b)/2.
    p3   = (log10(100.0/k)/d)
    p4   = (log10(eng/k)/d)
    return Amp*((eng/100.0)**p2)*(10**((p1*d*log((exp(p4)+exp(-p4))/2.)) \
                -(p1*d*log((exp(p3)+exp(-p3))/2.))))




def cutoffpl(energy, PhoIndex, HighECut, norm):
    '''
    cutoffpl(energy, PhoIndex, HighECut, norm)

    Parameters:
    -----------
    energy: float, energy to integrate over. 
    PhoIndex: float, power-law photon index. Same as alpha in band func,
                sbpl, etc. Essentailly the low-energy power-law index. 
                This value will have opposite sign as the RMFIT 'copl'.
    HighECut: float, e-folding energy of exponential rolloff (in keV).
                Same as E0 in band func. 
    norm: float, model amplitude. This is normalized to 1 keV where 
            RMFIT's is normalized to 100 keV. 

    XSPEC Cutoff Power-Law.
    https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node160.html

    
    Notes:
    ------
    This function is normalized to 1 keV and can not be changed. 

    This function is different from RMFIT's 'copl'. 
    This function has an opposite sign on the power law photon index
    and the function is normalized to 1 keV instead of 
    100 keV (like RMFIT's).
    To get RMFIT equivalent parameters:
    1.  Compute new Amplitude: 
            Amp(at 100 keV) = Amp*((1./100.)^PhoIndex)
    2.  Change sign on the power law photon index. 
    You should only do this for publishing the values. If you change the 
    parameters and then use this function later to compute flux, you will
    get the wrong flux. 

    norm * (energy**(-PhoIndex)) * exp(-energy/HighECut)

    eng = energy
    Amp = norm
    idx = PhoIndex
    HEC = HighECut

    Amp * (eng**(-idx)) * exp(-eng/HEC)

    '''
    eng = energy
    Amp = norm
    idx = PhoIndex
    HEC = HighECut
    return Amp * (eng**(-idx)) * exp(-eng/HEC)



def bbody(energy, kT, norm):
    '''
    bbody(energy, kT, norm)

    Parameters:
    ----------
    energy: float, energy to integrate over. 
    kT:     float, blackbody temperature parameter value. 
                Has units of keV. 
    norm:   float, normalization (aka amplitude) parameter value.  

    Notes:
    ------
    XSPEC Blackbody function.
    https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node137.html

    When integrating the blackbody function, you will run into the  
    "math domain error" 
    when kT is small. 

    For a GRB with z=0, kT < 14.1 keV would need a maximum limit.
    This is because 709.78*kT = 709.78*14.1 keV > 10,000 keV (10 MeV).
    We typically only integrate up to 10 MeV. 
    However, no GRBs will have z=0, this is just a limit. 

    emax = (709.7*kT)*(1+z)


    In order to avoid the math domain error, you should limit the 
    integration of the bbody up to 709.78 keV for every increase of 
    1 keV for kT. For example:
    kT     max energy
    --     ----------
    1       709.78
    2       1419.56
    3       2129.34
    10      7097.79
    The maximum allowed energy is 709.78 * kT.
    To be conservative, use 709.7 or even lower. 
    * Remember that energies will be redshifted (energy/(1+z)) when 
        calculating flux, fluence, and eiso. This max energy allowed is 
        the energy AFTER redshifting, not before. To account for this, 
        find what energy would be the max and give 709.7*kT in the 
        rest-frame of the burst. emax = (709.7*kT)*(1+z)


    The bbody will no longer be contributing to the total flux long 
    before you hit the maximum energy. 
    If this limit becomes a problem, try reducing the function you are 
    integrating to one where bbody is removed (i.e., grbm instead of 
        grbm+bbody or grbm+lpow instead of grbm+bbody+lpow). Since the 
        bbody is not contributing at those higher energies, this is 
        acceptable. 

    '''
    eng     = energy
    kT      = kT
    Amp     = norm
    return Amp * ( ((eng**2)*8.0525)/((kT**4)*(exp(eng/kT)-1)))
    #if eng <= (709.666 * kT): # to avoid exp overflow error
    #    return N * (((eng**2)*(8.0525)) / ((kT**4) * (exp(eng/kT)-1)))
    #else:
    #    return 0



def powerlaw(energy, PhoIndex, norm):
    '''
    powerlaw(energy, PhoIndex, norm)
    
    Parameters:
    ----------
    energy: float, energy element from an array of energies over which 
            the function is integrated over. 
    index:  float, power-law index parameter value.  
    norm:   float, normalization (aka amplitude) parameter value.  

    Notes:
    -------
    XSPEC's powerlaw
    https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node211.html
    This function is normalized to 1 keV instead of 100 keV. 
    It is different from RMFIT's by this and the sign change on the 
    power-law index. We wrote a function 'lpow' to be exactly the same 
    as RMFIT's. You can use that within XSPEC too and is also provided 
    as a function in this file. 
    '''
    eng    = energy
    Amp    = norm
    idx    = PhoIndex
    return Amp * (eng**(-idx))



def lpow(energy, plIndex, norm):
    '''
    lpow(energy, plIndex, norm)

    Parameters:
    ----------
    energy: float, energy element from an array of energies over which 
            the function is integrated over. 
    plIndex:  float, power-law index parameter value.  
    norm:   float, normalization (aka amplitude) parameter value.  


    Notes:
    ------
    This function is normalized to 100 keV where the 'powerlaw' is 
    normalized to 1 keV. 
    This is our user-defined power-law function. 
    It is the same as RMFIT's.
    XSPEC's 'powerlaw' normalizes the function to 1 keV instead 
    of 100 keV. 
    This version is the same as RMFIT's. 
    '''
    eng    = energy
    Amp    = norm
    idx    = plIndex
    return Amp * ((eng/100.)**idx)

