from __future__ import division
import numpy as np
from numpy import exp, log, log10
from math import atanh

__all__ = ["lpow", "copl", "band", "sbpl", "bbody"]


def band(energy, alpha, beta, epeak, norm):
    """

    """
    eng = energy
    a = alpha
    Amp = norm
    b = beta
    epk = epeak
    enorm = 100.0  # Don't have the option to change this!

    condLO = eng < ((a-b) * epk)/(2.0+a)
    condUP = eng >= ((a-b) * epk)/(2.0+a)

    bandLO = lambda x: Amp * pow((x/enorm), a) * exp(-(2.+a)*x/epk)
    bandUP = lambda x: Amp * pow((x/enorm), b) * exp(b-a) * \
                        pow(((a-b)*epk)/(enorm*(2.+a)), a-b)

    result = np.piecewise(eng, [condLO, condUP], [bandLO, bandUP])
    return result



def bbody(energy, kT, norm):
    '''

    '''
    eng = energy
    Amp = norm
    temp = kT
    return Amp * pow(eng, 2) * pow(exp(eng/temp)-1, -1)


def lpow(energy, index, norm):
    '''
    LPOW(energy, index, norm)

    RMFIT power-law.

    Should be same as XSPEC's PL model.
    
    eqn = A*(x/Epiv)**index
    '''
    eng = energy
    a = index
    Amp = norm
    enorm = 100.0 # don't have the option to changed this.
    result = Amp * pow(eng/enorm, a)
    return result


def copl(energy, index, epeak, norm):
    '''
    COPL(energy, index, epeak, norm)

    RMFIT's copl model.
    '''
    eng = energy
    a = index
    epk = epeak
    Amp = norm
    enorm = 100.0  # don't have the option to changed this.
    return Amp * exp(-eng*(2.+a)/epk) * pow(eng/enorm, index)



def sbpl(energy, alpha, beta, ebreak, norm):
    '''
    SBPL(energy, alpha, beta, ebreak, norm)

    RMFIT's SBPL model. Should be same as XSPEC's, but XSPEC's 
    has more options for the characteristic energy term.


    Conversion between epeak and ebreak:
    epeak = (ebreak) * (10**(0.3 * atanh((a+b+4)/(a-b))))
    ebreak = (epeak) / (10**(0.3 * atanh((a+b+4)/(a-b))))
           where a and b are alpha and beta

    '''

    eng  = energy
    a    = alpha
    b    = beta
    d    = 0.3 # break scale, don't have the option to change this.
    k    = ebreak
    Amp  = norm
    enorm = enorm # don't have the option to change this.

    p1   = (b-a)/2.
    p2   = (a+b)/2.
    p3   = log10(enorm/k)/d
    p4   = log10(eng/k)/d

    # pow(10, d) is same as 10**d
    part1 = p1 * d * log((exp(p4) + exp(-p4))/2.)
    part2 = p1 * d * log((exp(p3) + exp(-p3))/2.)
    
    result = Amp * pow(eng/enorm, p2) * pow(10, part1-part2)
    return result
