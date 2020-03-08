# ["lpow", "powerlaw", "cutoffpl", "grbm", "sbpl", "bbody"]
# 
#   grbm(energy, alpha, beta, tem, norm):
#   sbpl(energy, alpha, beta, ebreak, norm):
#   cutoffpl(energy, PhoIndex, HighECut, norm):
#   bbody(energy, kT, norm):
#   powerlaw(energy, PhoIndex, norm):
#   lpow(energy, plIndex, norm):

from __future__ import division
import warnings


def get_parameter_names(name):
    # pass one component at a time.
    if name == 'grbm':
        'alpha, beta, tem, norm'
    alpha, beta, ebreak, norm
    PhoIndex, HighECut, norm
    kT, norm
    PhoIndex, norm
    plIndex, norm



def check_parameters_type(modelname, parameters):
    if isinstance(parameters, dict):
        # dict type
        pass
    elif isinstance(parameters, list):
        # list of parameters.
        msg = "Warning: parameters in list must be in order. "
              "If using an additive model, pass a nested list of parameters."
        warnings.warn(msg, UserWarning)
        # Check for nested lists. All will be True if you do. 
        are_nested = [isinstance(sub, list) for sub in parameters]

        # `any` throws a False if you forget [] around a set of parameters.
        # use `all`
        if len(modelname.split('+')) > 1 and all(are_nested):
            # everything is fine
            print('a')
        elif len(modelname.split('+')) > 1 and not all(are_nested):
            msg = '''Make sure you are passing nested lists for additive models. 
            E.g., [[-1.25, -2.56, 500.0, 0.012], [51.8, 0.23]]'''
            raise Exception(msg)
        else:
            # you must have a single model component.
            pass


        #raise Exception("'parameters' must be of dictionary type.")


class Parameters(object):
    def __init__(self, modelname):
        self.modelname = modelname

    def _grbm(self):
