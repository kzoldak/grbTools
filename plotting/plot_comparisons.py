"""
This file contains the Eiso and parameter comparisons that are shown in 
    my dissertation. These are the plots with the diagonal line and the 
    x and y axes show identical variables, but different results of those 
    variables. 

"""

from __future__ import division

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from Zoldak.math.tools import root_mean_square as rms

# from .plotting import journal()


__all__ = ['get_parameter_label', 'plot_diagonal', 'plot_stats', 
            'plot_comparison', 'get_axes_limits']


def get_parameter_label(parameter):
    """
    Options:
    'eiso', 'epeak', 'epeakRest', 'alpha', 'beta', 'norm', 'flux', 'fluence'
    """
    labels = {
             'eiso': r'$E_{iso}$ (erg)', 
             'epeak': r'$E_{pk}$ (keV)',
             'epeakRest': r'$E^*_{pk}$ (keV)',
             'alpha': r'$\alpha$', 
             'beta': r'$\beta$', 
             'norm': r'$A$ (ph cm$^{-2}$ s$^{-1}$ keV$^{-1}$)',
             'flux': r'flux (erg cm$^{-2}$ s$^{-1}$)',
             'fluence': r'fluence (erg cm$^{-2}$)'
            }
    return labels[parameter]


def plot_diagonal(axLims, ax=None, **pltKwargs):
    """
    Leave this separate, so it doesn't get plotted more than once. 
    """      
    if ax is None:
        ax = plt.gca()

    xline = np.linspace(axLims[0], axLims[1], 100)
    yline = 1*xline+0 
    if pltKwargs:
        PLT = ax.plot(xline, yline, **pltKwargs)
    else:
        PLT = ax.plot(xline, yline, color='black', ls='-', lw=2, alpha=0.35)
    return PLT


def plot_stats(n, delta, sigma, xloc, yloc, ax=None):
    """
    Plot of statistics to overlap on your comparison plot. 
    Leave this separate. 

    Parameters:
    -----------
    n : int, 
        number of data points.
    delta : float, 
        average of the differences. Must log values before taking differences. 
        deltas = log10(eiso_2) - log10(eiso_1) 
        delta = mean(deltas), and the eiso_1 and eiso_2 are two arrays of 
        eiso energies for identical events. 
    sigma : float, 
        root mean square of the differences: rms(deltas)
    xloc : float, 
        Location on plot along x-axis for text. 0 to 1. 
    yloc : float, 
        Location on plot along y-axis for text. 0 to 1. 

    """
    if ax is None:
        ax = plt.gca()
    
    PLT = ax.text(0.05, 0.93, 
                  r'$n$ \ \ \ \ \ $\sigma$ \ \ \ \ \ \ \ $\Delta$',
                  fontsize=13, fontweight='bold', color ='black', 
                  transform=ax.transAxes) 

    ax.text(xloc, yloc, 
                  r'%.3f \ \ %.3f'%(sigma, delta), 
                  fontsize=12, color='black', 
                  transform=ax.transAxes)

    ax.text(0.04, yloc, r'%i:'%(n), 
            fontsize=12, fontweight='bold', color ='black', 
            transform=ax.transAxes)  
    return PLT


def get_axes_limits(x, y, axBuffer=0.15):
    axLims = (np.min([x.min(), y.min()]), 
              np.max([x.max(), y.max()]))
    if axBuffer:
        buff = (axLims[1]-axLims[0])*axBuffer
        axLims = (axLims[0]-buff, axLims[1]+buff)
    return axLims


# ********************************************************************
# ********************************************************************

def plot_comparison(x, y, xerr=None, yerr=None, parameter=None, xaxlabel='', 
    yaxlabel='', axLims=None, axBuffer=0.15, ax=None, **pltKwgs):
    """
    Plot for comparing Eiso and parameters, thos in my dissertation. 
    The ones with diagonals and have the same x- and y-axes labels. 


    Parameters:
    -----------
    x : array,
        x-axis values. Should already be logged, if needed. 
    y : array, 
        y-axis values. Should already be logged, if needed. 
    xerr : array, 
        error on x. Should already be logged, if needed. 
        Should be ** confidence intervals **
        E.g.:
        xerr = [df['eiso_ci_lo'].apply(np.log10), 
                df['eiso_ci_up'].apply(np.log10)]
    yerr : array, 
        error on x. See xerr.
    parameter : str, optional. 
        Parameter name. Options are:  
        'eiso', 'epeak', 'epeakRest', 'alpha', 'beta', 'norm', 'flux', 'fluence'
        If you do not wish to pass your own axlabel, use this attribute and 
          it will make one for you. Don't use this with axlabel. 
    xaxlabel : str, optional. 
        label for axes. Default is None. Don't use this with parameter. 
        It is recommended you use parameter, so that you get the right label. 
    yaxlabel : str, optional.
        See xaxlabel.
    axLims : list or tuple, optional. 
        x- and y- axes limits. (low, up). Default is None.
    axBuffer : float, optional. Range from 0 to 1, 0.15 is a good choice. 
        Provides an axis buffer of a fraction of the range. 0.15 is 15% of the 
        range of the data's min and max limits (for both axes). 
    ax : matplotlib axes object. Optional, default is None.
    **pltKwgs : dict, optional.
        Dictionary of plotting arguments for the data points. 
        Default is none. 


    Example:
    ---------
    x = df1['eiso'].apply(np.log10)
    y = df2['eiso'].apply(np.log10)
    
    xerr = [df1['eiso_ci_lo'].apply(np.log10), 
            df1['eiso_ci_up'].apply(np.log10)]
    
    yerr = [df2['eiso_ci_lo'].apply(np.log10), 
            df2['eiso_ci_up'].apply(np.log10)]
    
    axlabel = r'$E_{iso}$ (erg)'
    colors = ['blue',]
    labels=['Band vs SBPL',]
    
    plt.clf()
    plot_differences(x=x, y=x, xerr=xerr, yerr=yerr, axlabel=axlabel, 
                         axLims=None, ax=None)
    plt.show()
    
    
    """
    plt.rcParams.update( {'figure.subplot.bottom': 0.14,
                          'figure.subplot.left': 0.15,
                          'figure.subplot.right': 0.96,
                          'figure.subplot.top': 0.96,
                          'figure.subplot.wspace': 0.2,
                          'figure.subplot.hspace': 0.2,
                          })
    x = np.array(x)
    y = np.array(y)

    # if parameter and axlabel:
    #     raise AttributeError("Can't use both 'parameter' and 'axlabel'. Choose one.")
    # else:
    #     pass

    if parameter:
        parlabel = get_parameter_label(parameter)
    
    if ax is None:
        ax = plt.gca()
        
    if axLims is None:
        axLims = get_axes_limits(x=x, y=y, axBuffer=axBuffer)
    if not pltKwgs:
        pltKwgs = dict(marker='o', color='white', mec='blue', 
                       ms=6, lw=0, mew=0.9, alpha=1) # data points
    
    PLT = ax.plot(x, y, **pltKwgs)
    if xerr is not None:
        ax.hlines(x, xerr[0], xerr[1], color=pltKwgs['mec'])
    ax.set_xlim(*axLims)
    ax.set_ylim(*axLims)

    ax.set_xlabel('%s   (%s)'%(parlabel, xaxlabel), fontsize=12)

    labelpad = 0.0

    if ('flux' in parameter) or ('norm' in parameter):
        labelpad = -0.01

    ax.set_ylabel('%s   (%s)'%(parlabel, yaxlabel), fontsize=12, labelpad=labelpad)
    #ax.subplots_adjust(left=0.15, right=0.96, top=0.96, bottom=0.14)
    #ax.legend(loc=0, ) #fontsize=11, labelspacing=0.7, handletextpad=0, frameon=False)
    return PLT

# ********************************************************************
# ********************************************************************





# ax.set_xlabel('%s   (%s)'%(parlab, labels[0]), fontsize=12)
# if ('flux' in parlab) or ('norm' in parlab):
#     ax.set_ylabel('%s   (%s)'%(parlab, labels[1]), fontsize=12, labelpad=-0.01)
# else:
#     ax.set_ylabel('%s   (%s)'%(parlab, labels[1]), fontsize=12)
# return PLT





#     axLims = (np.min([x.min(), y.min()]), 
#               np.max([x.max(), y.max()]))
#     if axBuffer:
#         buff = (axLims[1]-axLims[0])*axBuffer
#         axLims = (axLims[0]-axbuff, axLims[1]+axbuff)
# Plot the diagonal line where x and y values are equal. 
# PLT = plot_diagonal(axLims, ax=ax)

# # Plot a line at equal values of x and y
# xlinedata = np.linspace(axLims[0], axLims[1], 10)
# ylinedata = 1.0*xlinedata+0
# PLT = ax.plot(xlinedata, ylinedata, 'k-', lw=2, alpha=0.8)
# for x,y,color,label in zip(xaxes,yaxes,colors,labels):


