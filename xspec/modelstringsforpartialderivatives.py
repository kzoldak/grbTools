'''
This file holds all the XSPEC model partial derivatives functions to 
be passed to sympy so that they can be computed for each model parameter. 
These are NOT the pre-computed partial derivatives. 

** It is recommended to use the pre-computed partial derivatives 
    since re-computing them every time is computationally expensive. 
    These pre-computed partial derivatives can be found in 
    modelpartialderivatives.py 
    of the same directory.
    **


These partial derivatives need to be calculated for flux, fluence, and 
eiso uncertainty, as well as epeak uncertainty when a model returns 
ebreak or tem or highECut instead of epeak. 

The model functions are set up with the same names as the XSPEC models:
'lpow', 'bbody', 'cutoffpl', 'grbm', 'sbpl'

**'lpow' is actually our user-defined powerlaw model** 
** 'sbpl' is our user-defined model**


Additive models:
  There is an '_' in place of '+' for additive models. 
  e.g., 'grbm_bbody_lpow' is the 'grbm+bbody+lpow' model. 


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




EXAMPLES:
---------
grbm_partials = epeak_partials('grbm')




def grbm(dpar, *parnames):
    pars        = list(parnames) + ['energy']
'''
from __future__ import division
from collections import OrderedDict
from grbTools.Math.partialderivatives import PartialDerivatives


def powerlaw(dpar=None):
    pars    = ['PhoIndex__1', 'norm__2', 'energy']
    eqn     = 'energy * (norm__2 * (energy**(-PhoIndex__1)))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut

def lpow(dpar=None):
    pars    = ['plIndex__1', 'norm__2', 'energy']
    eqn     = 'energy * (norm__2 * ((energy/100.)**plIndex__1))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut

def bbody(dpar=None):
    pars    = ['kT__1', 'norm__2', 'energy']
    eqn     = ('energy * (norm__2 * (((energy**2)*(8.0525))/'
               '((kT__1**4) * (exp(energy/kT__1)-1))))')
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut

def cutoffpl(dpar=None):
    '''
    The eqn should be:
    'energy * (norm__3 * ((energy/100.)**PhoIndex__1) * (exp(-energy/HighECut__2)))'

    This is the derivative of:
    'energy * (norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2)))'

    '''
    pars    = ['PhoIndex__1', 'HighECut__2', 'norm__3', 'energy']
    eqn     = 'energy * (norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2)))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut

def grbm(dpar=None):
    """
    The resulting dictionary will have two results in a list for each key.
    [lower, upper]
    """
    pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 'energy']
    lower = ('norm__4 * energy * ((energy/100.)**alpha__1)*'
             '(exp(-(energy/tem__3)))')
    upper = ('norm__4 * energy * ((((alpha__1-beta__2)*'
            'tem__3)/100.)**(alpha__1-beta__2)) * '
            '(exp(beta__2-alpha__1))*((energy/100.)**beta__2)')
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
        out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
        return out_L, out_U
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
            out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
            allOut[str(dpar)] = [out_L, out_U]
        return allOut


def sbpl(dpar=None):
    pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 'energy']
    eqn = 'energy * (norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp( log(energy/ebreak__3, 10)/0.3 ) + exp(-( log(energy/ebreak__3, 10)/0.3 )))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp(log(100./ebreak__3, 10)/0.3) + exp(-( log(100./ebreak__3, 10)/0.3 )))/2.)))))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut

def grbm_lpow(dpar=None):
    pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 
            'plIndex__5', 'norm__6',
            'energy']
    lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * ((energy/100.)**plIndex__5)))'
    upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + ((norm__6 * ((energy/100.)**plIndex__5))))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
        out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
        return out_L, out_U
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
            out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
            allOut[str(dpar)] = [out_L, out_U]
        return allOut


def grbm_powerlaw(dpar=None):
    pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 
            'PhoIndex__5', 'norm__6',
            'energy']
    lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (energy**(-PhoIndex__5))))'
    upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (energy**(-PhoIndex__5))))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
        out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
        return out_L, out_U
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
            out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
            allOut[str(dpar)] = [out_L, out_U]
        return allOut


def grbm_bbody(dpar=None):
    pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 
            'kT__5','norm__6',
            'energy']     
    lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))))'
    upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
        out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
        return out_L, out_U
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
            out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
            allOut[str(dpar)] = [out_L, out_U]
        return allOut

def grbm_bbody_lpow(dpar=None):
    pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 
            'kT__5','norm__6',
            'plIndex__7', 'norm__8',
            'energy']
    lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * ((energy/100.)**plIndex__7)))'
    upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + ((norm__8 * ((energy/100.)**plIndex__7))))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
        out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
        return out_L, out_U
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
            out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
            allOut[str(dpar)] = [out_L, out_U]
        return allOut


def grbm_bbody_powerlaw(dpar=None):
    pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 
            'kT__5','norm__6',
            'PhoIndex__7', 'norm__8',
            'energy']
    lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * (energy**(-PhoIndex__7))))'
    upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * (energy**(-PhoIndex__7))))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
        out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
        return out_L, out_U
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out_L   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
            out_U   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
            allOut[str(dpar)] = [out_L, out_U]
        return allOut


def cutoffpl_lpow(dpar=None):
    '''
    see cutoffpl.__doc__
    '''
    cutoffpl_lpow.__doc__ = cutoffpl.__doc__
    pars    = ['PhoIndex__1', 'HighECut__2', 'norm__3', 
               'plIndex__4', 'norm__5',
               'energy']
    eqn     = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * ((energy/100.)**plIndex__4)))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut

def cutoffpl_powerlaw(dpar=None):
    '''
    see cutoffpl.__doc__
    '''
    cutoffpl_lpow.__doc__ = cutoffpl.__doc__
    pars    = ['PhoIndex__1', 'HighECut__2', 'norm__3', 
               'PhoIndex__4', 'norm__5',
               'energy']
    eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (energy**(-PhoIndex__4))))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut


def cutoffpl_bbody(dpar=None):
    '''
    see cutoffpl.__doc__
    '''
    pars    = ['PhoIndex__1', 'HighECut__2', 'norm__3', 
               'kT__4', 'norm__5',
               'energy']
    eqn     = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (((energy**2)*(8.0525)) / ((kT__4**4) * (exp(energy/kT__4)-1)))))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut


def cutoffpl_bbody_lpow(dpar=None):
    '''
    see cutoffpl.__doc__
    '''
    pars    = ['PhoIndex__1', 'HighECut__2', 'norm__3', 
               'kT__4', 'norm__5',
               'plIndex__6', 'norm__7',
               'energy']
    eqn     = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (((energy**2)*(8.0525)) / ((kT__4**4) * (exp(energy/kT__4)-1)))) + (norm__7 * ((energy/100.)**plIndex__6)))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut

def cutoffpl_bbody_powerlaw(dpar=None):
    '''
    see cutoffpl.__doc__
    '''
    pars    = ['PhoIndex__1', 'HighECut__2', 'norm__3', 
               'kT__4', 'norm__5',
               'PhoIndex__6', 'norm__7',
               'energy']
    eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (((energy**2)*(8.0525)) / ((kT__4**4) * (exp(energy/kT__4)-1)))) + (norm__7 * (energy**(-PhoIndex__6))))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut

def sbpl_bbody(dpar=None):
    pars    = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 
               'kT__5', 'norm__6',
               'energy']
    eqn     = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut

def sbpl_lpow(dpar=None):
    pars    = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 
               'plIndex__5', 'norm__6',
               'energy']
    eqn     = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * ((energy/100.)**plIndex__5)))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut

def sbpl_powerlaw(dpar=None):
    pars    = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 
               'PhoIndex__5', 'norm__6',
               'energy']
    eqn     = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (energy**(-PhoIndex__5))))'    
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut


def sbpl_bbody_lpow(dpar=None):
    pars    = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 
               'kT__5', 'norm__6',
               'plIndex__7', 'norm__8',
               'energy']
    eqn     = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * ((energy/100.)**plIndex__7)))'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut


def sbpl_bbody_powerlaw(dpar=None):
    pars    = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 
               'kT__5', 'norm__6',
               'PhoIndex__7', 'norm__8',
               'energy']
    eqn = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * (energy**(-PhoIndex__7))))'    
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut



def epeak_partials(model, dpar=None):
    """

    Parameters:
    -----------
    model: str, model name. 
            Options are: 'grbm', 'cutoffpl', 'sbpl'
            THESE ARE XSPEC ONLY!
    dpar: str or None, default is None.
            Pass None and get a dictionary of all the partial derivatives.
            Pass a single string and get only the partial derivative
              wrt that parameter. 

    Notes:
    ------
    The 'sbpl' model uses brkscale = 0.3. If you don't want to fix it 
    at this value, update this function and add brkscale to the parameters.
    We currently don't allow the function to use this parameter, so for 
    now it doesn't matter. 

    """
    pars = None
    eqn = None
    if model == 'grbm':
        pars = 'alpha__1 tem__3'.split(' ')
        eqn  = '(alpha__1 + 2.) * tem__3'
        # pars = 'alpha tem'.split(' ')
        # eqn  = '(alpha+2.) * tem'
    elif model == 'cutoffpl':
        pars = 'PhoIndex__1 HighECut__2'.split(' ')
        eqn  = '(PhoIndex__1 + 2.) * HighECut__2'
        # pars = 'PhoIndex HighECut'.split(' ')
        # eqn  = '(PhoIndex+2.) * HighECut'
    elif model == 'sbpl':
        pars = 'alpha__1 beta__2 ebreak__3'.split(' ')
        eqn  = ('ebreak__3 * (10**(0.3 * (atanh('
                '(alpha__1 + beta__2 + 4)/(alpha__1 - beta__2)))))')
        # pars = 'alpha beta ebreak brkscale'.split(' ')
        # eqn  = 'ebreak*(10**(brkscale * (atanh((alpha + beta + 4)/(alpha-beta)))))'
    else:
        raise Exception, "Model not recognized."
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        for dpar in pars:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut



def ebreak_partials(dpar=None):
    """
    Maybe move this to the general directory. 
    This is from Kaneko et al 2006 appendix.

    """
    pars    = ['alpha', 'beta', 'epeak', 'energy']
    eqn     = '((alpha-beta)/2.)*(epeak/(2.+alpha))+12.5'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
        return out
    else:
        # Returns all partial derivatives in an OrderedDict.
        allOut = OrderedDict()
        # Don't need a partial derivative wrt 'energy'.
        for dpar in pars[:-1]:
            out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut


"""

Results of epeak_from_ebreak(dpar)
wrt each parameter:

partial deriv wrt:  alpha
10**(brkscale*atanh((alpha + beta + 4)/(alpha - beta)))*brkscale*ebreak*(1/(alpha - beta) - (alpha + beta + 4)/(alpha - beta)**2)*log(10)/(1 - (alpha + beta + 4)**2/(alpha - beta)**2)

partial deriv wrt:  beta
10**(brkscale*atanh((alpha + beta + 4)/(alpha - beta)))*brkscale*ebreak*(1/(alpha - beta) + (alpha + beta + 4)/(alpha - beta)**2)*log(10)/(1 - (alpha + beta + 4)**2/(alpha - beta)**2)

partial deriv wrt:  ebreak
10**(brkscale*atanh((alpha + beta + 4)/(alpha - beta)))

partial deriv wrt:  brkscale
10**(brkscale*atanh((alpha + beta + 4)/(alpha - beta)))*ebreak*log(10)*atanh((alpha + beta + 4)/(alpha - beta))

"""


# def epeak_from_ebreak(dpar=None):
#     """
#     epeak_from_ebreak(dpar=None)

#     If dpar=None, runs all partial derivatives in default.
#     If you pass a parameter to dpar, it will comute only that partial 
#       derivative. 
#     dpar options:  'alpha' 'beta' 'ebreak' and 'brkscale'

#     When using 'sbpl' model.
#     We set brkscale to 0.3 when modeling, so we really don't need to 
#       take the derivative of it. But we leave it in this form for now.  
#     """
#     pars = 'alpha beta ebreak brkscale'.split(' ')
#     eqn  = 'ebreak*(10**(brkscale * (atanh((alpha + beta + 4)/(alpha-beta)))))'
    
#     if dpar is not None:
#         out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
#         return out
#     else:
#         allOut = OrderedDict()
#         for dpar in pars:
#             out  = str(PartialDerivatives(eqn, str(dpar), 0, *pars))
#             allOut[str(dpar)] = out
#         return allOut
