"""
Only the most general math related tools should go here. 
Those tools specific to certain tasks associated with my dissertation should 
go within the subdirectores. 
"""
import numpy as np

__all__ = ['weighted_average', 'percentage_increase', 'percentage_decrease', 
    'percentage_difference', 'percentage_error', 'residuals', 'root_mean_square', 
    'sum_of_squares', 'to_log']

def weighted_average(values, ndatas):
    """
    Takes the weighted average of a list of values, given the number of 
    ndata points. 

    Parameters:
    -----------
    values:  list, 
        of values that need to be weighted by the
                number of data points in their respective samples. 
    ndatas:  list, 
        of the number of data points in the samples. 
    
    Returns:
    --------
    Returns the individual weights in a list and then the sum of the 
    weights (which is the weighted average).

    We pull the N into the calculation of each ith value since this 
    gives us the individual weights. However, this equation is 
    typically shown with N pulled outside. 

    WA = (1/N)*SUM(n_i * v_i)  
    where n_i is the number of data points in each ith iteration and
    v_i is the value associated with the ith iteration. 
    If this were weighted averages of sigma:
    WA = (1/N) * SUM(n_i * sigma_i)

    n = [13, 32, 35]
    sigma = [0.015, 0.104, 0.096]
    then N = 13+32+35 = 80
     
    
    """
    N = np.sum(ndatas)
    ##          sigma*(n/N)
    weights = [(v*(n/N)) for n,v in zip(ndatas, values)]
    return weights, np.sum(weights)


def percentage_increase(original_value, new_value):
    '''
    Increase in a value by percentage. 

    Parameters:
    -----------
    original_value : int or float. Original number. 
    new_value : int or float. New number. 

    Returns:
    --------
    Increase of new value over the original value, as a percentage. 
    
    Equation:
    ---------
    ((new_value - original_value)/original_value) * 100.0

    See Also:
    --------
    percentage_decrease
    percentage_error
    percentage_difference
    
    '''
    return ((new_value - original_value)/original_value) * 100.0


def percentage_decrease(original_value, new_value):
    '''
    Decrease in a value by percentage. 

    Parameters:
    -----------
    original_value : int or float. Original number. 
    new_value : int or float. New number. 

    Returns:
    --------
    Float, dncrease of new value relative to the original value, as a percentage. 
    
    Equation:
    ---------
    ((original_value - new_value)/original_value) * 100.0

    See Also:
    --------
    percentage_increase
    percentage_error
    percentage_difference

    
    '''
    return ((original_value - new_value)/original_value) * 100.0



def percentage_error(actual, estimated):
    '''
    Percentage error equation.

    Parameters:
    -----------
    actual : float, actual value. 
    estimated : float, estimate of the actual value. 

    Returns:
    --------
    Float

    Equation:
    ---------
    difference = (abs(estimated - actual)/(actual)) * 100.0

    Notes:
    -----
    Actual would be a theoretical value and estiamted would be 
    the experimental value. If you don't have a *true* actual value, as in 
    know or highly supported by theory, use percentage_difference instead.

    See Also:
    --------
    percentage_increase
    percentage_decrease
    percentage_difference
    
    '''
    difference = (abs(estimated - actual)/(actual)) * 100.0
    return difference



def percentage_difference(a, b):
    '''
    Parameters:
    -----------
    a : list or array, of values. Indices between a and b must match up.
    b : list or array, of values. Indices between a and b must match up.

    Math:
    -----
    abs((abs(b - a) / (0.5*(b + a)) ) * 100.0)

    **Percentage differences should always be positive.

    Returns:
    --------
    An array of percentage differences between the elements within 
    array a and array b. 
    
    Notes:
    ------
    Percentage Difference calculation between two arrays, a and b.
    
    Each of the ith elements in the arrays (or lists, if you pass lists) 
     must match up. E.g., a[0] will be compared against b[0], and so on.
     If you pass lists, they will be converted to arrays and returned as 
     an array of percentage differences. Arrays are easier to work with 
     for mathematical calculations. 

    See Also:
    --------
    percentage_increase
    percentage_decrease
    percentage_error
    '''
    a = np.asarray(a)
    b = np.asarray(b)
    return abs((abs(b - a) / (0.5*(b + a)) ) * 100.0)


    
def residuals(ydata, ymodel):
    '''
    Parameters:
    -----------
    ydata : list or array. 
        True y-axis values to compare against the model determined values. 
    ymodel : list or array. 
        Model determined y-axis values. These lie perfectly on the model curve. 

    Returns:
    --------
    List,
        of residuals, which are the differences between the true data and the 
        model: data - model
        These are unweighted residuals. 
        The sigma residuals are weighted: (data-model)/sigma
        This function does not do this. 
    
    Equation:
    ---------
    resids = [(i-j) for i,j in zip(ydata-ymodel)]
    '''
    resids = [(i-j) for i,j in zip(ydata-ymodel)]
    return resids
    

def sum_of_squares(values):
    """
    Sum of the squared values in a list. 

    Parameters:
    -----------
    values:  list or array,
         of values to square and then sum.

    Returns:
    --------
    List, of the sum of the squared values.

    Equation:
    ---------
    sum( [i**2 for i in values] )

    """
    return sum( [i**2 for i in values] )


def root_mean_square(values):
    """
    Root mean square of the differences, where 'values' is a list of 
    differences. 


    Parameters:
    -----------
    values:  list,
         of differences between two sets of values to take the rms of.

    Equation:
    ---------
    np.sqrt( sum_of_squares(values) / len(values) )
    where 
    sum_of_squares(values) equation is:
    sum( [i**2 for i in values] )
    """
    #return np.sqrt(sum_of_squares(values)/len(values))
    return (sum_of_squares(values)/len(values))**0.5







def _check_error_type(x, xerr):
    """
    Checks error type; confidence intervals or margins of error. 
    
    Returns
    --------
    'ci' if your errors are confidence intervals.
    'moe' if your errors are margins of error. 
    
    Function raises an excpetion if your errors are bad. 
    - 
    
    """
    x = np.asarray(x)
    xerr = np.asarray(xerr)
    if xerr.ndim != 2:
        raise Exception('xerr should have 2 dimensions, [ErrLo, errUp]')
    xerrLo = xerr[0]
    xerrUp = xerr[1]
    if (any(x > xerrLo) or any(xerrUp > x)) and (all(x > xerrLo) and all(xerrUp > x)):
        # You have confidence intervals
        print('You have Confidence Intervals')
        return('ci')
    elif all(abs(xerrLo) < abs(x)) and all(abs(xerrUp) < abs(x)):
        # all errors lower than the values indicates margins of error. 
        print('You have Margins of Error')
        return('moe')
    else:
        print('You have problems with your errors. ')
        raise Exception('You have problems with your errors.')
        

        
def to_log(x, xerr, which='both', errTypeReturn='moe'):
    """
    Take linear data and uncertainties and transform them to log
    values. 
    
    Parameters
    ----------
    x : array of floats
        measurements of which to take logarithms
        
    xerr : array of floats
        uncertainties on x. Can pass errors with either 1 or 2 dimensions. 
        Note that 1d errors will be assumed to be Margins of Error and symmetric 
        when linear. Symmetry is not preserved when logging, so output will be 
        slightly asymmetric. 
        
    which : str, 'lower', 'upper', 'both', or 'average'
        Which uncertainty to return. Default is 'both'.
        
        Note that when converting to/from
        linear and logarithmic spaces, errorbar symmetry is not
        preserved. You can no log margins of error. You must first 
        find the confidence intervals, log the values and the confidence 
        intervals, then take the differences to get the logged margins 
        of error. 
        
    errTypeReturn : str, 'moe' or 'ci'.
        Default is 'moe'. 
        'moe' returns logged margins of error.
        'ci' returns logged confidence intervals. 
        Note that which='average' can not be used with errTypeReturn='ci' 
        since you can not have averaged confidence intervals. 
        

    Returns
    -------
    logx : array of floats
    logxLo : array of floats
    logxUp : array of floats
    
    Notes:
    ------
    Note that 1d errors (xerr) will be assumed to be Margins of Error and symmetric 
        when linear. Symmetry is not preserved when logging, so output will be 
        slightly asymmetric. 
    
    This function uses another function to check whether your xerr are Confidence Intervals 
        or Margin of Error. 
    '_check_error_type(x, xerr=xerr)' returns a string of either 'ci' or 'moe'.
    
    """
    x = np.asarray(x)
    xerr = np.asarray(xerr)
    
    if xerr.ndim == 1:
        # If xerr is only 1 dimension, then you MUST have Margin of Error (MOE). 
        # Can't log MOE, convert to Confidence Intervals (CI) before logging.
        logxLo, logx, logxUp = np.log10([x-xerr, x, x+xerr])
        
    elif xerr.ndim == 2:
        # If xerr is 2 dimensions, you could have MOE or CI. 
        errtype = _check_error_type(x, xerr=xerr)  # check error type. 
        if errtype == 'moe':
            # If you have margins of error.
            moeLo,moeUp = xerr
            logxLo, logx, logxUp = np.log10([x-moeLo, x, x+moeUp])
        else:
            # If you have confidence intervals.
            ciLo,ciUp = xerr
            logxLo, logx, logxUp = np.log10([ciLo, x, ciUp])
    else:
        raise Exception('xerr must be 1 or 2 dimensions.')
    
    # logxLo and logxUp are currently in Conf Intv form. Convert to Marg of Err here, if desired.
    
    if errTypeReturn == 'moe':
        # logxLo and logxUp are currently in Conf Intv form. Convert to Marg of Err here, if desired.
        logxLo = logx-logxLo
        logxUp = logxUp-logx     
        if which == 'average':
            return logx, 0.5*(logxLo+logxUp)
        elif which == 'lower':
            return logx, logxLo
        elif which == 'upper':
            return logx, logxUp
        else:
            return logx, logxLo, logxUp
    if errTypeReturn == 'ci':
        if which == 'average':
            raise Exception("Can't return Confidence Intervals that are averaged.")
        elif which == 'lower':
            return logx, logxLo
        elif which == 'upper':
            return logx, logxUp
        else:
            return logx, logxLo, logxUp

