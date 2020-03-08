"""
This file contains all of the propagation of uncertainty tools for the 
XSPEC modeling functions. 

"""

from __future__ import division
from collections import OrderedDict
import json

__all__ = ['partial_derivative', 'xspec_str_functions', 
            'calc_partial_derivatives', 'calc_partial_derivative', 
            'model_partial_derivatives']


def partial_derivative(func, dpar, parameters, pretty_print=False):
    """
    Parameters
    ----------
    func :   str 
        The function's equation as a string. 
                **SEE NOTES (below) ON LOGS**
    dpar :   str 
        Parameter you wish to take the derivative of func with respect to. 
    parameters : list of str
        A list of the parameter names in str format. These parameter names 
        should match those used in func. See Equation and parameters (below).
    pretty_print : True, False, optional
        Default is False. 
        True turns on pretty printing so that equations are output to the 
        screen in a readable format. 
        True is recommended within Jupyter notebooks. 
    
    Returns
    -------
    str
    A string format of the partial derivative wrt a single parameter, 
        as defined by dpar.

    Notes
    -----
    !!!! WARNING ON SYMPY LOGS !!!!
    Sympy requires input of ln for log and log for log10, however, it 
        outputs log for ln and log10 for log base 10, which is what you 
        are used to using in Python. So the output can be used immediatly 
        within Python. 
        For every one of the functions in this file, we replaced all log 
        with ln and all log10 with log.
    
    

    Examples
    --------

    """
    import sympy
    from sympy import Function, Symbol
    from sympy import integrate, diff, exp, log, ln, sqrt, lambdify
    from sympy import init_printing
    from sympy.parsing.sympy_parser import parse_expr
    from sympy import mpmath
    #from sympy.mpmath import atanh
    from sympy import atanh
    
    init_printing(pretty_print=pretty_print) # True turns it on.
    
    for par in parameters:
        # This includes 'energy'.
        locals()[par] = Symbol('%s'%par, real=True)
        
    return diff(eval(func), eval(dpar), method='quad')



def xspec_str_functions(modelname):
    # Any modelname having 'grbm' in it will have to run the 
    #   partial_derivative function twice.
    modelOptions = ['powerlaw', 'lpow', 'bbody', 
                    'grbm', 'sbpl', 'cutoffpl', 
                    'grbm+lpow', 'sbpl+lpow', 'cutoffpl+lpow', 
                    'grbm+powerlaw', 'sbpl+powerlaw', 'cutoffpl+powerlaw',
                    'grbm+bbody', 'sbpl+bbody', 'cutoffpl+bbody', 
                    'grbm+bbody+lpow', 'sbpl+bbody+lpow', 
                    'cutoffpl+bbody+lpow', 
                    'grbm+bbody+powerlaw', 'sbpl+bbody+powerlaw', 
                    'cutoffpl+bbody+powerlaw', 
                    'bbody+lpow', 'bbody+powerlaw']
    
    if modelname == 'powerlaw':
        eqn = 'energy * (norm__2 * (energy**(-PhoIndex__1)))'
        pars = ['PhoIndex__1', 'norm__2']
        
    elif modelname == 'lpow':
        eqn = 'energy * (norm__2 * ((energy/100.)**plIndex__1))'
        pars = ['plIndex__1', 'norm__2']
        
    elif modelname == 'bbody':
        eqn = ('energy * (norm__2 * (((energy**2)*(8.0525))/'
               '((kT__1**4) * (exp(energy/kT__1)-1))))')
        pars = ['kT__1', 'norm__2']
        
    elif modelname == 'cutoffpl':
        eqn = ('energy * (norm__3 * ((energy)**-PhoIndex__1) '
               '* (exp(-energy/HighECut__2)))')
        pars = ['PhoIndex__1', 'HighECut__2', 'norm__3']
        
    elif modelname == 'grbm':
        pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4']
        lower = ('norm__4 * energy * ((energy/100.)**alpha__1)*'
                 '(exp(-(energy/tem__3)))')
        upper = ('norm__4 * energy * ((((alpha__1-beta__2)*'
                'tem__3)/100.)**(alpha__1-beta__2)) * '
                '(exp(beta__2-alpha__1))*((energy/100.)**beta__2)')
        eqn = [lower, upper]
        
    elif modelname == 'sbpl':
        pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4']
        eqn = ('energy * (norm__4*((energy/100.)**((alpha__1+beta__2)/2.))'
               '* (10**((((beta__2 - alpha__1)/2.) * 0.3 * '
               'ln((exp(log(energy/ebreak__3, 10)/0.3 ) + '
               'exp(-( log(energy/ebreak__3, 10)/0.3 )))/2.)) - '
               '(((beta__2 - alpha__1)/2.) * 0.3 * '
               'ln((exp(log(100./ebreak__3, 10)/0.3) + '
               'exp(-( log(100./ebreak__3, 10)/0.3 )))/2.)))))')
    
    elif modelname == 'grbm+lpow':
        pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 'plIndex__5', 'norm__6']
        lower = ('energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) '
                 '+ (norm__6 * ((energy/100.)**plIndex__5)))')
        upper = ('energy * ((norm__4 *((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2))'
                 ' * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + '
                 '((norm__6*((energy/100.)**plIndex__5))))')
        eqn = [lower, upper]
    
    elif modelname == 'grbm+powerlaw':
        pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 'PhoIndex__5', 'norm__6']
        lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (energy**(-PhoIndex__5))))'
        upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (energy**(-PhoIndex__5))))'
        eqn = [lower, upper]
        
    elif modelname == 'grbm+bbody':
        pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 'kT__5','norm__6']    
        lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))))'
        upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))))'
        eqn = [lower, upper]
        
    elif modelname == 'grbm+bbody+lpow':
        pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 'kT__5','norm__6','plIndex__7', 'norm__8']
        lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * ((energy/100.)**plIndex__7)))'
        upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + ((norm__8 * ((energy/100.)**plIndex__7))))'
        eqn = [lower, upper]
    
    elif modelname == 'grbm+bbody+powerlaw':
        pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4','kT__5','norm__6','PhoIndex__7', 'norm__8']
        lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * (energy**(-PhoIndex__7))))'
        upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * (energy**(-PhoIndex__7))))'
        eqn = [lower, upper]
    
    elif modelname == 'cutoffpl+lpow':
        pars = ['PhoIndex__1', 'HighECut__2', 'norm__3', 'plIndex__4', 'norm__5']
        eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * ((energy/100.)**plIndex__4)))'
    
    elif modelname == 'cutoffpl+powerlaw':
        pars = ['PhoIndex__1', 'HighECut__2', 'norm__3', 'PhoIndex__4', 'norm__5']
        eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (energy**(-PhoIndex__4))))'
    
    elif modelname == 'cutoffpl+bbody':
        pars = ['PhoIndex__1', 'HighECut__2', 'norm__3', 'kT__4', 'norm__5']
        eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (((energy**2)*(8.0525)) / ((kT__4**4) * (exp(energy/kT__4)-1)))))'
    
    elif modelname == 'cutoffpl+bbody+lpow':
        pars = ['PhoIndex__1', 'HighECut__2', 'norm__3', 'kT__4', 'norm__5',
                'plIndex__6', 'norm__7']
        eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (((energy**2)*(8.0525)) / ((kT__4**4) * (exp(energy/kT__4)-1)))) + (norm__7 * ((energy/100.)**plIndex__6)))'
    
    elif modelname == 'cutoffpl+bbody+powerlaw':
        pars = ['PhoIndex__1', 'HighECut__2', 'norm__3', 'kT__4', 'norm__5',
                'PhoIndex__6', 'norm__7']
        eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (((energy**2)*(8.0525)) / ((kT__4**4) * (exp(energy/kT__4)-1)))) + (norm__7 * (energy**(-PhoIndex__6))))'
    
    elif modelname == 'sbpl+bbody':
        pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 'kT__5', 'norm__6']
        eqn = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))))'
    
    elif modelname == 'sbpl+lpow':
        pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 'plIndex__5', 'norm__6']
        eqn = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * ((energy/100.)**plIndex__5)))'
    
    elif modelname == 'sbpl+powerlaw':
        pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 'PhoIndex__5', 'norm__6']
        eqn = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (energy**(-PhoIndex__5))))'    
    
    elif modelname == 'sbpl+bbody+lpow':
        pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 'kT__5', 'norm__6',
               'plIndex__7', 'norm__8']
        eqn = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * ((energy/100.)**plIndex__7)))'
    
    elif modelname == 'sbpl+bbody+powerlaw':
        pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 
                'kT__5', 'norm__6',
                'PhoIndex__7', 'norm__8']
        eqn = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * (energy**(-PhoIndex__7))))'    
    
    elif modelname == 'bbody+lpow':
        eqn = ('energy * ((norm__2 * (((energy**2)*(8.0525))/((kT__1**4)*(exp(energy/kT__1)-1)))) '
                '+ (norm__4 * ((energy/100.)**plIndex__3)))')   
        pars = ['kT__1', 'norm__2', 'plIndex__3', 'norm__4']
        
    elif modelname == 'bbody+powerlaw':
        eqn = ('energy * ((norm__2 * (((energy**2)*(8.0525))/((kT__1**4)*(exp(energy/kT__1)-1)))) '
                '+ (norm__4 * (energy**(-PhoIndex__3))))')        
        pars = ['kT__1', 'norm__2', 'PhoIndex__3', 'norm__4']
    
    else:
        msg = ("'modelname' not recognized. Must be one of following: "
               "%r"%(modelOptions))
        raise Exception(msg)
    
    if 'energy' in pars:
        raise Exception('remove energy from pars')
    
    return eqn, pars

        

def calc_partial_derivatives(modelname):
    """
    This function computes ALL partial derivatives. 
    """

    equation,parameters = xspec_str_functions(modelname)
    
    partials = OrderedDict()
    for dpar in parameters:
        if isinstance(equation, list):
            # for models with 'grbm' in them, there is a lower and upper part.
            out = []
            for eqn in equation:
                pderiv = str(partial_derivative(func=eqn, 
                                                dpar=str(dpar), 
                                                parameters=parameters+['energy'], 
                                                pretty_print=False))
                out.append(pderiv)
        else:
            # all other models
            out = str(partial_derivative(func=equation, 
                                         dpar=str(dpar), 
                                         parameters=parameters+['energy'], 
                                         pretty_print=False))
        # Store each one in OrderedDict.
        partials[str(dpar)] = out
    return partials


def calc_partial_derivative(modelname, dpar):
    """
    This function computes a single partial derivative. 
    """

    equation,parameters = xspec_str_functions(modelname)
    
    if isinstance(equation, list):
        # for models with 'grbm' in them, there is a lower and upper part.
        out = []
        for eqn in equation:
            pderiv = str(partial_derivative(func=eqn, 
                                            dpar=str(dpar), 
                                            parameters=parameters+['energy'], 
                                            pretty_print=False))
            out.append(pderiv)
    else:
        # all other models
        out = str(partial_derivative(func=equation, 
                                     dpar=str(dpar), 
                                     parameters=parameters+['energy'], 
                                     pretty_print=False))
    return out


def model_partial_derivatives(modelname, dpar=None):
    f = ('/Users/KimiZ/Python/My_Modules/grbTools/xspec/'
         'xspec_functions_partial_derivatives.json')
    partials = json.load(open(f, 'r'), object_pairs_hook=OrderedDict)
    
    if dpar is not None:
        return partials[modelname][dpar]
    return partials[modelname]


def epeak_partials(modelname, dpar=None):
    """
    Parameters
    ----------
    modelname : str
        Model name in string format. 
        Funtion recognizes 'grbm', 'cutoffpl', or 'sbpl' in the model's 
        name and will calculate epeak unc accordingly. 
    dpar : str, None, optional
        Parameter name, in string format, that you wish to compute the 
        partial derivative of. 
        None is default and returns a dict of all of them. 

    Notes
    -----
    We will always want the Epeak of the PRIMARY model component, and 
        Band, Sbpl, or Cutoffpl will ALWAYS be the primary component. 
        Therefore the *__1, *__2, etc in the parameter names work, 
        regardless of whether we have a single model or three additive 
        models combined into one. 

    Examples
    --------
    epeak_partials(modelname='grbm', dpar=None)
    epeak_partials(modelname='grbm', dpar='alpha__1')
    epeak_partials(modelname='cutoffpl', dpar=None)
    epeak_partials(modelname='sbpl', dpar=None)
    """
    if 'grbm' in modelname:
        pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4']
        eqn  = '(alpha__1 + 2.) * tem__3'

    elif 'cutoffpl' in modelname:
        pars = ['PhoIndex__1', 'HighECut__2', 'norm__3']
        eqn  = '(PhoIndex__1 + 2.) * HighECut__2'

    elif 'sbpl' in modelname:
        pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4']
        eqn  = ('ebreak__3 * (10**(0.3 * (atanh('
                '(alpha__1 + beta__2 + 4.)/(alpha__1 - beta__2)))))')
    else:
        raise Exception("'modelname' not recognized.")
    
    if dpar is not None:
        out = str(partial_derivative(func=eqn, 
                                     dpar=str(dpar), 
                                     parameters=pars, 
                                     pretty_print=False))
        return out
    else:
        allOut = OrderedDict()
        for dpar in pars:
            out = str(partial_derivative(func=eqn, 
                                     dpar=str(dpar), 
                                     parameters=pars, 
                                     pretty_print=False))
            allOut[str(dpar)] = out
        return allOut


def epeak_partial_derivatives(modelname, dpar=None):
    f = ('/Users/KimiZ/Python/My_Modules/grbTools/xspec/'
         'xspec_epeak_partial_derivatives.json')
    partials = json.load(open(f, 'r'), object_pairs_hook=OrderedDict)
    if dpar is not None:
        return partials[modelname][dpar]
    return partials[modelname]
       


#  THIS IS NOT READY TO USE YET, DON'T INCLUDE IN THE __all__ list!!
def ebreak_partials(dpar=None):
    """
    Maybe move this to the general directory. 
    This is from Kaneko et al 2006 appendix.

    """
    pars    = ['alpha', 'beta', 'epeak']
    eqn     = '((alpha-beta)/2.)*(epeak/(2.+alpha))+12.5'
    if dpar is not None:
        # Returns ONLY the one partial derivative. 
        out  = str(calc_partial_derivative(eqn, str(dpar), 0, *pars))
        return out
    else:
        allOut = OrderedDict()
        for dpar in pars[:-1]:
            out  = str(calc_partial_derivative(eqn, str(dpar), 0, *pars))
            allOut[str(dpar)] = out
        return allOut 
