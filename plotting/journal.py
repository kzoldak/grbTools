"""
Only the most basic plotting functions should go here. 
If you need to, you can place others here but they should be moved to 
another subdirectory that is more appropriate for its task. 

"""

from __future__ import division
import matplotlib.pyplot as plt


__all__ = ['journal',]


def journal():
    '''
    Setup for consistent, journal worthy matplotlib figures. 
  
    Notes:
    ------
    Takes no arguments, simply run: 
    journal()
    Performs:  plt.rcParams.update(params)    
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



