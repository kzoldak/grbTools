"""
The tools included in this file are specific to comparing values of Eiso 
energies or spectral modeling parameters. Most functions are statistical in nature. 



"""
import numpy as np
import pandas as pd

from grbTools.math.tools import root_mean_square as rms


__all__ = ['statistical_uncertainty', 'weighted_average_uncertainty', 
    'uncertainty_ratios', 'quadrature_uncertainties', 
    'compare_param_uncertainty_overlap', 'parameter_pretty_print']

# *********************************************************
# **************    STATISTICAL TOOLS   *******************
# *********************************************************

def statistical_uncertainty(dataframe, parameter, suffix='', logit=False):
    """
    dataframe : pandas dataframe. Finds the values and confidence intervals 
                    within the dataframe. 

    parameter : str, string name of the parameter. 
                e.g., 'alpha', 'eiso', ...

    suffix : str, default is ''. If the parameter name and its errors has a 
                suffix of '_x' or '_y' from joining dataframes, place that 
                in the suffix. When we merge dataframes, our '_x' and '_y' 
                are after the '_moe_lo', so adding '_x' at the end of the 
                parameter string wont work.

    logit : True or False, default is False. 
            logit = True will apply log10 to the data and confidence intervals.
    
    Prints:
    -------
    Float. The sample's statistical uncertainty, rounded to 3 decimal places.
    
    Returns:
    -------
    Float. Same as prints, but does not round.
    
    Math:
    -----
    Averages the upper and lower confidence intervals (in their margin of 
    error form) for each GRB in the sample and then takes the average of the 
    all of those to give the *sample* average. 
    
    
    Example:
    --------
    for par,logit in zip(['eiso', 'alpha', 'beta', 'epeak', 'norm'],
                         [True, False, False, True, True]):
        statistical_uncertainty(dataframe=df1, parameter=par, logit=logit)
    
    """
    # only a single dataframe
    d = dataframe
    
    if logit == True:
        values = d[parameter+'%s'%suffix].apply(np.log10)
        ci_lo = d[parameter+'_ci_lo'+'%s'%suffix].apply(np.log10)
        ci_up = d[parameter+'_ci_up'+'%s'%suffix].apply(np.log10)
    else:
        values = d[parameter+'%s'%suffix]
        ci_lo = d[parameter+'_ci_lo'+'%s'%suffix]
        ci_up = d[parameter+'_ci_up'+'%s'%suffix] 
    # Three ways to do this. 
    #a = np.mean([np.mean([i,j]) for i,j in zip(values-ci_lo, ci_up-values)])
    #b = np.mean([values-ci_lo, ci_up-values])
    #c = (pd.DataFrame([values-ci_lo, ci_up-values]).mean()).mean()
    #[print('%.3f'%i) for i in [a,b,c]];
    a = np.mean([np.mean([i,j]) for i,j in zip(values-ci_lo, ci_up-values)])
    print('%.3f'%a)
    return a


def weighted_average_uncertainty(sigmas, nGrbs):
    """
    Parameters:
    -----------
    sigmas : list of array. The systematical uncertainties of 
                of each study.
    nGrbs : list or array. The number of grbs in each of the studies that found the above sigmas.
            For each ith sigma, the ith nGrb must match. 
            For example, if the ith sigma is 0.014 and 28 GRBs were included in the study to find 
            that value of 0.014, then nGrb=28.
            
    Returns
    -------
    weights, weighted uncertainties, sum of the weighted uncertainties.

    Notes:
    ------
    This is the same as the weighted arithmetic mean formula, except instead of
    measurements of x we have systematical uncertainties (the sigmas). 
    Our weights are given by: nGrbs[i]/sum(nGrbs), 
    the number of GRBs in the ith place over the sum of all the GRBs included 
    in the weighted average calculation. 
    Our weights sum to 1.
            
    """
    nGrbs = np.asarray(nGrbs)
    sigmas = np.asarray(sigmas)
    N = sum(nGrbs)
    weights = []
    weighted_sigmas = []
    for n,sig in zip(nGrbs, sigmas):
        weight = n/N
        weighted_sigma = weight*sig
        weights.append(weight)
        weighted_sigmas.append(weighted_sigma)
    print(weights)
    print('')
    print(weighted_sigmas)
    print('')
    print('\nWeights')
    print(['%.3f'%i for i in weights])
    print('\nWeighted Sigmas')
    print(['%.3f'%i for i in weighted_sigmas])
    print('\nWeighted Average Sigma')
    print('%.3f'%sum(weighted_sigmas))
    print('\n')
    return weights, weighted_sigmas, sum(weighted_sigmas)
 

def uncertainty_ratios(sys, stats):
    """
    Parameters:
    -----------
    sys : float. Systematical uncertainty.
    stats : list or array of the sample statistical uncertainties.
            There should only be 2 since 2 data sets were 
            involved in calculaing sys.
    
    Returns:
    --------
    Prints and returns the fractional ratios of sys/stat for each stat 
        in stats.

    """
    stats = np.asarray(stats)
    ratios = sys/stats
    print('\nsys/stat ratios')
    print(['%.3f'%i for i in ratios])
    return ratios


def quadrature_uncertainties(sys, stats):
    """
    Parameters:
    -----------
    sys : float. Systematical uncertainty.
    stats : list or array of the sample statistical uncertainties.
            There should only be 2 since 2 data sets were 
            involved in calculaing sys.
    
    Returns:
    --------
    Two parts. 
    Prints and returns the statistical and systematical uncertainties 
    added in quadrature.
        (sys**2 + stat**2)**0.5
    """
    stats = np.asarray(stats)
    quad_results = np.asarray([(i**2 + sys**2)**0.5 for i in stats])
    print('\n sqrt(sys^2 + stat^2) -- uncertainties added in quadrature.')
    print(['%.3f'%i for i in quad_results])
    return quad_results


def compare_param_uncertainty_overlap(valuesA, valuesB, errorsA, errorsB):
    """
    Compare two sets of values and their 1-sigma margins of error for 
    overlap.  If value returned is between 0 and 1, they overlap. 

    *** 
        THIS FUNCTION GOES ALONG WITH 
        plot_param_uncertainty_overlap(valuesA, valuesB, errorsA, errorsB, deltaValues)
        IN  plotting.py
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
    deltaValues = abs(valuesA - valuesB)/(errorsA + errorsB)  
    return deltaValues


# *********************************************************
# **************     PARAMETER TOOLS    *******************
# *********************************************************

def parameter_pretty_print(parameter):
    """
    Parameters:
    -----------
    parameter : str, parameter name. 
                options (see Notes).

    Notes:
    ------
    Given the parameter name, this function will return the proper string 
    for printing consistent values of that parameter. For example, if you 
    want alpha, beta, or z values, you typically want them rounded to 3 
    decimal places. This function would return '%.3f' for any of those three
    parameter names. 
    For epeak or epeakRest, this function returns '%.0f'. Using '%i' will 
    improperly round the values. 
    For 'norm', 'eiso', 'flux', 'fluence', this function will return '%.3E'.

    This function takes the following parameter names:
    'alpha', 'beta', 'z'
    'epeak', 'epeakRest'
    'norm', 'eiso', 'flux', 'fluence'

    Example:
    --------

    for par in ['alpha', 'beta', 'epeak', 'epeakRest', 'norm', 'eiso', 'z']:
        form = pretty_print(par)
        printstr = '%s & %s & %s \\\\'%(par, form, form ) 
        print(printstr%(df[par].mean(), df[par].median()))
    
    
    """
    if parameter in ['alpha', 'beta', 'z']:
        return '%.3f'
    elif parameter in ['epeak', 'epeakRest']:
        return '%.0f'
    elif parameter in ['norm', 'eiso', 'flux', 'fluence']:
        return '%.3E'





    