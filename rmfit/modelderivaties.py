from __future__ import division
# from collections import OrderedDict
import os, sys, json
# import pandas as pd
# import numpy
from numpy import log, log10, exp, sqrt
# from scipy import integrate

from Zoldak.Math.partialderivatives import PartialDerivatives



old_stdout = sys.stdout  # SAVE ORIGINAL PRINT TO SCREEN SETTINGS.
# DISABLE PRINTING
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
# RESTORE PRINTING
def enablePrint():
    sys.stdout = old_stdout


all_models = ['band', 'band+bbody', 'grbm+lpow', 'grbm+bbody+lpow', 'sbpl', 'sbpl+bbody', 'sbpl+lpow', 'sbpl+bbody+lpow', 'cutoffpl', 'cutoffpl+bbody', 'cutoffpl+lpow', 'cutoffpl+bbody+lpow']




def band(energy):
    '''
    THIS ONE WORKS FOR RMFIT. DO NOT EDIT!!!!
    
    '''
    pars        = ['alpha__1','beta__2','epeak__3','norm__4','energy'] 
    if energy < ((alpha__1 - beta__2) * (epeak__3/(2. + alpha__1))):
        lower = 'norm__4 * energy * ((energy/100.)**alpha__1)*(exp(-(energy/(epeak__3/(2. + alpha__1)))))'
        out   = str(PartialDerivatives(lower, str(dpar), 0, *pars))
        return eval(out)
    else:
        upper = 'norm__4 * energy * ((((alpha__1-beta__2)*(epeak__3/(2. + alpha__1)))/100.)**(alpha__1-beta__2)) * (exp(beta__2-alpha__1))*((energy/100.)**beta__2)'
        out   = str(PartialDerivatives(upper, str(dpar), 0, *pars))
        return eval(out) 