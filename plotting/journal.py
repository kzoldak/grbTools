from __future__ import division

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


__all__ = ['journal',]


def journal():
    '''

    Journal worthy setting for matplotlib plots. 
    Takes no arguments; simply run: 
    journal()


    Notes:
    ------

    Appears nearly square. I like this best.
    'figure.figsize': [3.1, 2.6]
    
    Use: plt.tight_layout(pad=0.1, w_pad=0.0, h_pad=0.0)
    plt.xlabel('$E_{iso}$ $(erg)$',labelpad=-1)  
    plt.ylabel('$E^*_{pk}$ $(keV)$',labelpad=-2)
    
    '''
    params = {'backend': 'pdf',
              'axes.labelsize':  10,
              'font.size':       10,
              'legend.fontsize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'xtick.direction': 'in',
              'ytick.direction': 'in',
              'text.usetex':     True,
              #'text.usetex':     False,
              'figure.figsize': [4,3], #[4,3], #[3.1, 2.6], # [7,6]
              'font.family': 'serif',
              'ytick.right': True}
    plt.rcParams.update(params)



# def journal():
#     '''
#     Journal worth setting for matplotlib plots. 
#     Takes no arguments; simply run: 
#     journal()
#     Notes:
#     ------
#     Appears nearly square. I like this best.
#     'figure.figsize': [3.1, 2.6]
    
#     Use: plt.tight_layout(pad=0.1, w_pad=0.0, h_pad=0.0)
#     plt.xlabel('$E_{iso}$ $(erg)$',labelpad=-1)  
#     plt.ylabel('$E^*_{pk}$ $(keV)$',labelpad=-2)
    
#     '''
#     params = {'backend': 'pdf',
#               'axes.labelsize':  10,
#               'font.size':       10,
#               'legend.fontsize': 8,
#               'xtick.labelsize': 8,
#               'ytick.labelsize': 8,
#               'ytick.right': True,
#               'xtick.direction': 'in',
#               'ytick.direction': 'in',
#               'text.usetex':     True,
#               #'text.usetex':     False,
#               'figure.figsize': [4,3], #[4,3], #[3.1, 2.6], # [7,6]
#               'font.family': 'serif',}
#     plt.rcParams.update(params)


    