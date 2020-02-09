#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 11:23:14 2020

@author: kimzoldak
"""



import numpy as np

from numpy import exp

import matplotlib.pyplot as plt
%matplotlib inline

import astropy
from astropy.modeling import Fittable1DModel, Parameter, fitting


from lumdist import Cosmology as Cosmo
cosmo = Cosmo()


def calc_reduced_chi_square(fit, x, y, yerr, N, n_free):
    '''
    fit (array) values for the fit
    x,y,yerr (arrays) data
    N total number of points
    n_free number of parameters we are fitting
    '''
    return 1.0/(N-n_free)*sum(((fit - y)/yerr)**2)


class FLCDM(Fittable1DModel):
    """
    All parameters go here.
    """
    H0 = Parameter(default=1.)
    Om = Parameter(default=1.)
    Ol = Parameter(default=1.)

    @staticmethod
    def evaluate(x, H0, Om, Ol):
        c = 2.99792458E+5
        return np.log10((c*(1.+x)/H0) * (((1.+x)**3.)*Om + Ol)**-0.5)

    @staticmethod
    def fit_deriv(x, H0, Om, Ol):
        """
        Logged DL.
        These are the derivatives of the DL func for the Lambda CDM model. 
        x - the z that is not part of the integral.
        """
        c = 2.99792458E+5
        #d_H0 = -0.434294/H0
        #d_H0 = np.ones(x.shape)
        d_H0 = np.ones(x.shape) * -0.434294/H0
        d_Om = -(0.217147*(x + 1)**3)/(Om*(x + 1)**3 + Ol)
        d_Ol = -0.217147/(Om*(x + 1)**3 + Ol)
        return [d_H0, d_Om, d_Ol]




# FlatLambdaCDM 
redshifts = np.linspace(0.1,8,100)
cosmoconstants = {'H0':67.8, 'Om': 0.308}
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='FlatLambdaCDM', 
                                 units='Mpc',
                                 **cosmoconstants)

# x2 : redshifts,  y2 : LumDist
x2 = redshifts
y2 = np.log10(DL)
np.random.seed(1) #seed
y2 += np.random.normal(0., 0.1, x2.shape) # 0.5
y2_err = np.ones(x2.shape)*0.1


plt.clf()
plt.errorbar(x2, y2, yerr=y2_err, fmt='.')
plt.show()



m = FLCDM(H0=75., Om=0.27, Ol=0.73)
fitter = fitting.LevMarLSQFitter()
f = fitter(model=m, x=x2, y=y2, weights=1.0/y2_err**2)


plt.errorbar(x2, y2, yerr=y2_err, fmt='.')
plt.plot(x2, f(x2))
plt.show()

#print(f)
#calc_reduced_chi_square(f(x2), x2, y2, y2_err, len(x2), 3)

#
#
#
#
#
#sine_model = SineNew(a=4.,b=2.,c=4.,d=0.)
#fitter = fitting.LevMarLSQFitter()
#sine_fit = fitter(sine_model, x2, y2, weights = 1.0/y2_err**2)
#
#
#
#plt.errorbar(x2, y2, yerr=y2_err, fmt='.')
#plt.plot(x2,sine_fit(x2))
#plt.show()
#
#print(sine_fit)
#
#calc_reduced_chi_square(sine_fit(x2), x2, y2, y2_err, len(x2), 3)
#

