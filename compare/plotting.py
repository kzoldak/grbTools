
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from grbTools.plotting.journal import journal
journal()


__all__ = ['plot_parameter_comparisons', 'plot_param_uncertainty_overlap']
        

def plot_parameter_comparisons(data1, data2):
    """
    Plotting function for side-by-side parameter values or energetics. 
    This is NOT the function for the dissertation plots with the diagonal 
    line. 


    Parameters:
    ------------
    data1 : three arrays of data:  values, lowerCIs, upperCIs.
    data2 : three arrays of data:  values, lowerCIs, upperCIs.
    
    Example:
    -------
    for param in ['alpha', 'beta', 'epeak', 'norm', 'eiso']:
        if param == 'eiso':
            dataA = list(zip(df1[param].apply(np.log10), 
                             df1[param+'_ci_lo'].apply(np.log10),
                             df1[param+'_ci_up'].apply(np.log10)))
            
            dataB = list(zip(df2[param].apply(np.log10), 
                             df2[param+'_ci_lo'].apply(np.log10),
                             df2[param+'_ci_up'].apply(np.log10)))
        else:
            dataA = list(zip(df1[param], 
                             df1[param+'_ci_lo'],
                             df1[param+'_ci_up']))
            
            dataB = list(zip(df2[param], 
                             df2[param+'_ci_lo'],
                             df2[param+'_ci_up']))
        plot_parameter_comparisons(dataA, dataB)  

    
    """
    data1 = np.asarray(data1)
    data2 = np.asarray(data2)
    v1 = data1[:,0]
    l1 = data1[:,1]
    u1 = data1[:,2]
    v2 = data2[:,0]
    l2 = data2[:,1]
    u2 = data2[:,2]

    def format_axes(fig):
        for i, ax in enumerate(fig.axes):
            #ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
            ax.tick_params(labelbottom=False, 
                           labeltop=False, 
                           labelleft=False, 
                           labelright=False, 
                           labelsize=None) 
            ax.set_xlim(-0.25, 0.5)
            ax.set_xlabel(None)
            ax.set_ylabel(None)
    fig,ax = plt.subplots(figsize=(16,2), nrows=1, ncols=29, 
                          sharex=False, 
                          sharey=False, squeeze=True)
    gs = GridSpec(1, 29, figure=fig, wspace=0, hspace=0)
    
    for i in range(0, 29):
        ax[i] = fig.add_subplot(gs[0, i:i+1])
        ax[i].plot(0, v1[i], 'o', color='red', lw=0, alpha=1)
        ax[i].vlines(0, l1[i], u1[i], color='red',  alpha=0.5)
        ax[i].plot(0.25, v2[i], 'o', color='blue', lw=0, alpha=1)
        ax[i].vlines(0.25, l2[i], u2[i], color='blue',  alpha=0.5)
    format_axes(fig)
    plt.show()



def plot_param_uncertainty_overlap(valuesA, valuesB, errorsA, errorsB, 
    deltaValues, burstID=None):
    """
    *** 
        THIS FUNCTION GOES ALONG WITH 
        compare_param_uncertainty_overlap(valuesA, valuesB, errorsA, errorsB)
        IN  tools.py
    ***


    Parameters:
    -----------
    valueA : array, of one set of values.
    valueB : array, of other set of values.
    errorsA : array, average errors on valueA (margins of error).  
                If asymmetrical, average them first.  
    errorsB : array, average errors on valueB (margins of error).  
                If asymmetrical, average them first.
    
    Notes:
    ------
    Make sure that the events being compared are IDENTICAL! 
        I.e., trigger or GRB number match for each event in both 
        valuesA and valuesB arrays. 
    
    This statistic renders a value of unity when the deviation between the
        two values (one from A and one from B arrays)
        exactly matches the sum of the 1-sigma errors. 
        A value < 1 indicates the values are within 1-sigma errors, and
        a value > 1 indicates that the values are NOT
        within 1-sigma errors of each other. 

    
    Returns:
    -------
    An array of values, explained in Notes. 

    """
    valuesA = np.asarray(valuesA)
    valuesB = np.asarray(valuesB)
    errorsA = np.asarray(errorsA)
    errorsB = np.asarray(errorsB)

    if burstID is not None:
        burstID = np.asarray([str(i+':  ') for i in burstID])
    else:
        burstID = np.asarray(['' for i in range(0, len(valuesA))])
        
    pltKwargsA = dict(fmt='ro', alpha=0.9, label='valuesA')
    pltKwargsB = dict(fmt='bo', alpha=0.9, label='valuesB')

    i = 0
    for vala,valb,erra,errb,bid in zip(valuesA, valuesB, \
                                        errorsA, errorsB, burstID):
        plt.clf()
        plt.figure(figsize=(4,5))
        plt.errorbar(x=0.5, y=vala, xerr=None, yerr=erra, **pltKwargsA)
        plt.errorbar(x=1, y=valb, xerr=None, yerr=errb, **pltKwargsB)
        plt.xlim(0, 1.5)
        if deltaValues[i] < 1.0:
            title = '** %.3f **'%deltaValues[i]
        else:
            title = '%.3f'%deltaValues[i]
        plt.title('%s%s'%(bid, title)) # this likely looks confusing. 
        plt.legend(loc=0)
        plt.show()
        i += 1

