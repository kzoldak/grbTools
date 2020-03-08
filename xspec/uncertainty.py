from __future__ import division
from collections import OrderedDict
import json

from partialderivatives import PartialDerivatives


deriv = PartialDerivatives()
deriv.calc_partial_derivatives(modelname='cutoffpl', 
                               dpar=None, 
                               epeakUnc=False)


def flux_uncertainty():
    # If covarmat available, do this.
    # else do this.
    pass

def epeak_uncertainty():
    # If covarmat available, do this.
    # else do this.
    pass