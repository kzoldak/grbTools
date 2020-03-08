"""
This file contains all of the propagation of uncertainty tools for the 
XSPEC modeling functions. 

"""
from __future__ import division
from collections import OrderedDict
import json
import os

import grbTools
# base location of module.
# grbTools.__path__[0]  # /Users/kimzoldak/Github/grbTools 


class PartialDerivatives(object):
    """
    This class propagates flux uncertainty from the parameter errors and the 
    covariance matrices from the XSPEC/PyXSPEC fits. 
    This provides errors on fluence and eiso as well. 
    This class also computes the epeak uncertainties. To invoke this feature, 
    use epeakunc=True. The default is False, so that flux uncertainties are 
    always computed otherwise. 

    We don't recommend runnung the 
    'calc_partial_derivative' 
    and 
    'calc_partial_derivatives' 
    functions over and over, as they employ the sympy package to compute the 
    partial derivatvies of the functions every time you want a flux or epeak 
    uncertainty. 
    Instead, use 
    ''
    and
    ''
    which read the precomputed partial derivatives from a JSON file. 


    Examples
    --------
    deriv = PartialDerivatives()

    deriv.modelnameOptions  # get modelname options.

    deriv.calc_partial_derivatives(modelname='sbpl', dpar=None, epeakUnc=False)

    deriv.calc_partial_derivatives(modelname='cutoffpl', dpar=None, epeakUnc=False)

    deriv.calc_partial_derivatives(modelname='grbm', dpar=None, epeakUnc=False)

    deriv.calc_partial_derivatives(modelname='grbm', 
                                   dpar=None, 
                                   epeakUnc=False)['alpha__1'][0]

    deriv.calc_partial_derivatives(modelname='grbm', 
                                 dpar='alpha__1', 
                                 epeakUnc=False)[0]

    deriv.calc_partial_derivatives(modelname='grbm', 
                                 dpar=None, 
                                 epeakUnc=False)['alpha__1'][1]

    deriv.calc_partial_derivatives(modelname='grbm', 
                                 dpar='alpha__1', 
                                 epeakUnc=False)[1]

    print(deriv.calc_partial_derivatives(modelname='grbm', 
                                   dpar=None, 
                                   epeakUnc=True))

    print(deriv.calc_partial_derivatives(modelname='grbm+bbody', 
                                   dpar=None, 
                                   epeakUnc=True))

    print(deriv.calc_partial_derivatives(modelname='grbm+bbody+lpow', 
                                   dpar=None, 
                                   epeakUnc=True))

    deriv.calc_partial_derivatives(modelname='sbpl', 
                                   dpar=None, 
                                   epeakUnc=True)

    deriv.calc_partial_derivatives(modelname='sbpl+bbody', 
                                   dpar=None, 
                                   epeakUnc=True)

    deriv.calc_partial_derivatives(modelname='sbpl+bbody+lpow', 
                                   dpar=None, 
                                   epeakUnc=True)

    deriv.calc_partial_derivatives(modelname='cutoffpl', 
                                   dpar=None, 
                                   epeakUnc=True)

    deriv.calc_partial_derivatives(modelname='cutoffpl+bbody', 
                                   dpar=None, 
                                   epeakUnc=True)

    deriv.calc_partial_derivatives(modelname='cutoffpl+bbody+lpow', 
                                   dpar=None, 
                                   epeakUnc=True)

    """
    # In future, consider inheriting from Directories and Files.
    def __init__(self, modelname=None, dpar=None, epeakUnc=False):
        self.modelnameOptions = ['powerlaw', 'lpow', 'bbody', 
                                'grbm', 'sbpl', 'cutoffpl', 
                                'grbm+lpow', 'sbpl+lpow', 'cutoffpl+lpow', 
                                'grbm+powerlaw', 'sbpl+powerlaw', 'cutoffpl+powerlaw',
                                'grbm+bbody', 'sbpl+bbody', 'cutoffpl+bbody', 
                                'grbm+bbody+lpow', 'sbpl+bbody+lpow', 
                                'cutoffpl+bbody+lpow', 
                                'grbm+bbody+powerlaw', 'sbpl+bbody+powerlaw', 
                                'cutoffpl+bbody+powerlaw', 
                                'bbody+lpow', 'bbody+powerlaw']
        if modelname:
            self.modelname = modelname
            if modelname not in self.modelnameOptions:
                msg = ("'modelname' not recognized. Must be one of following: "
                       "%r"%(self.modelnameOptions))
                raise Exception(msg)
        if dpar:
            self.dpar = dpar

        if epeakUnc:
            self.epeakUnc = epeakUnc
            options = ['grbm', 'sbpl', 'cutoffpl']
            if not any(map(lambda x: x in modelname, options)):
                msg2 = ("modelname must have 'grbm', 'sbpl' or 'cutoffpl' "
                        "in its name to get epeak uncertainty.")
                raise Exception(msg2)
        
    def partial_derivative(self, func, dpar, parameters, pretty_print=False):
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

        Imports the follwing functions only:
            atanh, exp, log, ln, sqrt


        Examples
        --------
        eqn = 'energy * (norm__2 * ((energy/100.)**plIndex__1))'
        pars = ['plIndex__1', 'norm__2']
        
        par = pars[0]

        out = partial_derivative(func = eqn, dpar = str(par), 
                                 parameters = pars+['energy'], 
                                 pretty_print = False)
        """
        # These imports must STAY local to this function. 
        from sympy import Symbol
        from sympy import diff
        from sympy import atanh, exp, log, ln, sqrt
        for par in parameters:      # This includes 'energy'.
            locals()[par] = Symbol('%s'%par, real=True)
        return diff(eval(func), eval(dpar), method='quad')


    def _xspec_str_model_functions(self):
        """
        Grab the parameters and model function in string format for 
        computing partial derivatives. 

        """
        # Any self.modelname having 'grbm' in it will have to run the 
        #   partial_derivative function twice.
        
        if self.modelname == 'powerlaw':
            eqn = 'energy * (norm__2 * (energy**(-PhoIndex__1)))'
            pars = ['PhoIndex__1', 'norm__2']
            
        elif self.modelname == 'lpow':
            eqn = 'energy * (norm__2 * ((energy/100.)**plIndex__1))'
            pars = ['plIndex__1', 'norm__2']
            
        elif self.modelname == 'bbody':
            eqn = ('energy * (norm__2 * (((energy**2)*(8.0525))/'
                   '((kT__1**4) * (exp(energy/kT__1)-1))))')
            pars = ['kT__1', 'norm__2']
            
        elif self.modelname == 'cutoffpl':
            eqn = ('energy * (norm__3 * ((energy)**-PhoIndex__1) '
                   '* (exp(-energy/HighECut__2)))')
            pars = ['PhoIndex__1', 'HighECut__2', 'norm__3']
            
        elif self.modelname == 'grbm':
            pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4']
            lower = ('norm__4 * energy * ((energy/100.)**alpha__1)*'
                     '(exp(-(energy/tem__3)))')
            upper = ('norm__4 * energy * ((((alpha__1-beta__2)*'
                    'tem__3)/100.)**(alpha__1-beta__2)) * '
                    '(exp(beta__2-alpha__1))*((energy/100.)**beta__2)')
            eqn = [lower, upper]
            
        elif self.modelname == 'sbpl':
            pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4']
            eqn = ('energy * (norm__4*((energy/100.)**((alpha__1+beta__2)/2.))'
                   '* (10**((((beta__2 - alpha__1)/2.) * 0.3 * '
                   'ln((exp(log(energy/ebreak__3, 10)/0.3 ) + '
                   'exp(-( log(energy/ebreak__3, 10)/0.3 )))/2.)) - '
                   '(((beta__2 - alpha__1)/2.) * 0.3 * '
                   'ln((exp(log(100./ebreak__3, 10)/0.3) + '
                   'exp(-( log(100./ebreak__3, 10)/0.3 )))/2.)))))')
        
        elif self.modelname == 'grbm+lpow':
            pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 'plIndex__5', 'norm__6']
            lower = ('energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) '
                     '+ (norm__6 * ((energy/100.)**plIndex__5)))')
            upper = ('energy * ((norm__4 *((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2))'
                     ' * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + '
                     '((norm__6*((energy/100.)**plIndex__5))))')
            eqn = [lower, upper]
        
        elif self.modelname == 'grbm+powerlaw':
            pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 'PhoIndex__5', 'norm__6']
            lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (energy**(-PhoIndex__5))))'
            upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (energy**(-PhoIndex__5))))'
            eqn = [lower, upper]
            
        elif self.modelname == 'grbm+bbody':
            pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 'kT__5','norm__6']    
            lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))))'
            upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))))'
            eqn = [lower, upper]
            
        elif self.modelname == 'grbm+bbody+lpow':
            pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4', 'kT__5','norm__6','plIndex__7', 'norm__8']
            lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * ((energy/100.)**plIndex__7)))'
            upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + ((norm__8 * ((energy/100.)**plIndex__7))))'
            eqn = [lower, upper]
        
        elif self.modelname == 'grbm+bbody+powerlaw':
            pars = ['alpha__1', 'beta__2', 'tem__3', 'norm__4','kT__5','norm__6','PhoIndex__7', 'norm__8']
            lower = 'energy * ((norm__4 * ((energy/100.)**alpha__1)*(exp(-(energy/tem__3)))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * (energy**(-PhoIndex__7))))'
            upper = 'energy * ((norm__4 * ((((alpha__1-beta__2)*tem__3)/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * (energy**(-PhoIndex__7))))'
            eqn = [lower, upper]
        
        elif self.modelname == 'cutoffpl+lpow':
            pars = ['PhoIndex__1', 'HighECut__2', 'norm__3', 'plIndex__4', 'norm__5']
            eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * ((energy/100.)**plIndex__4)))'
        
        elif self.modelname == 'cutoffpl+powerlaw':
            pars = ['PhoIndex__1', 'HighECut__2', 'norm__3', 'PhoIndex__4', 'norm__5']
            eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (energy**(-PhoIndex__4))))'
        
        elif self.modelname == 'cutoffpl+bbody':
            pars = ['PhoIndex__1', 'HighECut__2', 'norm__3', 'kT__4', 'norm__5']
            eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (((energy**2)*(8.0525)) / ((kT__4**4) * (exp(energy/kT__4)-1)))))'
        
        elif self.modelname == 'cutoffpl+bbody+lpow':
            pars = ['PhoIndex__1', 'HighECut__2', 'norm__3', 'kT__4', 'norm__5',
                    'plIndex__6', 'norm__7']
            eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (((energy**2)*(8.0525)) / ((kT__4**4) * (exp(energy/kT__4)-1)))) + (norm__7 * ((energy/100.)**plIndex__6)))'
        
        elif self.modelname == 'cutoffpl+bbody+powerlaw':
            pars = ['PhoIndex__1', 'HighECut__2', 'norm__3', 'kT__4', 'norm__5',
                    'PhoIndex__6', 'norm__7']
            eqn = 'energy * ((norm__3 * ((energy)**-PhoIndex__1) * (exp(-energy/HighECut__2))) + (norm__5 * (((energy**2)*(8.0525)) / ((kT__4**4) * (exp(energy/kT__4)-1)))) + (norm__7 * (energy**(-PhoIndex__6))))'
        
        elif self.modelname == 'sbpl+bbody':
            pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 'kT__5', 'norm__6']
            eqn = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))))'
        
        elif self.modelname == 'sbpl+lpow':
            pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 'plIndex__5', 'norm__6']
            eqn = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * ((energy/100.)**plIndex__5)))'
        
        elif self.modelname == 'sbpl+powerlaw':
            pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 'PhoIndex__5', 'norm__6']
            eqn = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (energy**(-PhoIndex__5))))'    
        
        elif self.modelname == 'sbpl+bbody+lpow':
            pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 'kT__5', 'norm__6',
                   'plIndex__7', 'norm__8']
            eqn = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * ((energy/100.)**plIndex__7)))'
        
        elif self.modelname == 'sbpl+bbody+powerlaw':
            pars = ['alpha__1', 'beta__2', 'ebreak__3', 'norm__4', 
                    'kT__5', 'norm__6',
                    'PhoIndex__7', 'norm__8']
            eqn = 'energy * ((norm__4 * ((energy/100.)**((alpha__1 + beta__2)/2.)) * (10.**((((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(energy/ebreak__3, 10)/0.3)) + exp(-(log(energy/ebreak__3, 10)/0.3)))/2.)) - (((beta__2 - alpha__1)/2.) * 0.3 * ln((exp((log(100./ebreak__3, 10)/0.3)) + exp(-(log(100./ebreak__3, 10)/0.3)))/2.))))) + (norm__6 * (((energy**2)*(8.0525)) / ((kT__5**4) * (exp(energy/kT__5)-1)))) + (norm__8 * (energy**(-PhoIndex__7))))'    
        
        elif self.modelname == 'bbody+lpow':
            eqn = ('energy * ((norm__2 * (((energy**2)*(8.0525))/((kT__1**4)*(exp(energy/kT__1)-1)))) '
                    '+ (norm__4 * ((energy/100.)**plIndex__3)))')   
            pars = ['kT__1', 'norm__2', 'plIndex__3', 'norm__4']
            
        elif self.modelname == 'bbody+powerlaw':
            eqn = ('energy * ((norm__2 * (((energy**2)*(8.0525))/((kT__1**4)*(exp(energy/kT__1)-1)))) '
                    '+ (norm__4 * (energy**(-PhoIndex__3))))')        
            pars = ['kT__1', 'norm__2', 'PhoIndex__3', 'norm__4']

        else:
            msg = ("'modelname' not recognized. Must be one of following: "
                   "%r"%(self.modelnameOptions))
            raise Exception(msg)
        if 'energy' in pars:
            raise Exception('remove energy from pars')
        return eqn, pars


    def _xspec_str_epeak_functions(self):
        """
        Grab the parameters and the epeak eqn in string format for 
        computing partial derivatives.  
        """
        pars = self._xspec_str_model_functions()[1]  # grab parameters only
        if 'grbm' in self.modelname:
            eqn  = '(alpha__1 + 2.) * tem__3'
        elif 'sbpl' in self.modelname:
            eqn = ('ebreak__3 * (10**(0.3 * (atanh('
                   '(alpha__1 + beta__2 + 4.)/(alpha__1 - beta__2)))))')
        elif 'cutoffpl' in self.modelname:
            eqn  = '(PhoIndex__1 + 2.) * HighECut__2'
        return eqn, pars


    def calc_partial_derivatives(self, modelname, dpar=None, epeakUnc=False):
        """
        This function computes ALL partial derivatives, but if dpar is 
        given, then only one of them is returned. We retired the option of 
        only computing a single partial derivative.

        Uses external function: 
        'partial_derivative'
        """
        
        PartialDerivatives.__init__(self, modelname, dpar, epeakUnc)
        
        if epeakUnc is True:
            # Equation here will always be a str. 
            eqn,pars = self._xspec_str_epeak_functions() # epeak fn
        else:
            # Equation could be a list of the model has 'grbm' in it.
            eqn,pars = self._xspec_str_model_functions() # model fn

        partials = OrderedDict()
        for par in pars:
            if isinstance(eqn, list):
                # Models with 'grbm' will be a list with 2 parts:
                #   [0] is lower Band condition and [1] is the upper. 
                eqns = eqn # dont need copy bc we aren't chaning it as we iter.
                out = []
                for eqn in eqns:
                    pderiv = self.partial_derivative(func = eqn, 
                                                     dpar = par, 
                                                     parameters = pars+['energy'], 
                                                     pretty_print = False)
                    out.append(str(pderiv))
            else:
                # all other models, only a single eqn
                pderiv = self.partial_derivative(func = eqn, 
                                                 dpar = par, 
                                                 parameters = pars+['energy'], 
                                                 pretty_print = False)
                out = str(pderiv)
            partials[par] = out # Store all Part Derivs in an OrderedDict.
        if dpar is not None:
            return partials[str(dpar)]  # only one partial deriv is returned.
        return partials  # all are retuned as an OrderedDict


    def get_partial_derivatives(self, modelname, dpar=None, epeakUnc=False):
        PartialDerivatives.__init__(self, modelname, dpar, epeakUnc)
        path = grbTools.__path__[0]
        if epeakUnc is True:
            f = os.path.join(path, 'xspec/files/xspec_epeak_partial_derivatives.json')
            # '/Users/KimiZ/Python/My_Modules/grbTools/xspec/'
            #  'xspec_epeak_partial_derivatives.json')
            partials = json.load(open(f, 'r'), object_pairs_hook=OrderedDict)
            if dpar is not None:
                return partials[modelname][dpar]
            return partials[modelname]
        else:
            f = os.path.join(path, 'xspec/files/xspec_functions_partial_derivatives.json')
            # f = ('/Users/KimiZ/Python/My_Modules/grbTools/xspec/'
            #      'xspec_functions_partial_derivatives.json')
            partials = json.load(open(f, 'r'), object_pairs_hook=OrderedDict)
            if dpar is not None:
                return partials[modelname][dpar]
            return partials[modelname]


# def model_partial_derivatives(modelname, dpar=None):
#     f = ('/Users/KimiZ/Python/My_Modules/grbTools/xspec/'
#          'xspec_functions_partial_derivatives.json')
#     partials = json.load(open(f, 'r'), object_pairs_hook=OrderedDict)
#     if dpar is not None:
#         return partials[modelname][dpar]
#     return partials[modelname]


# def epeak_partial_derivatives(modelname, dpar=None):
#     f = ('/Users/KimiZ/Python/My_Modules/grbTools/xspec/'
#          'xspec_epeak_partial_derivatives.json')
#     partials = json.load(open(f, 'r'), object_pairs_hook=OrderedDict)
#     if dpar is not None:
#         return partials[modelname][dpar]
#     return partials[modelname]



