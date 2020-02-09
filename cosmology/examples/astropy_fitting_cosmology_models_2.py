#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 11:07:41 2020

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


a : H0
m : Om
y : Ol


da
dm
dy


d/da((c*(1.+x)/a) * (((1.+z)**3.)*m + y)**-0.5) 


d_da = -(c*(x + 1))/(H0**2 * (Om*(z + 1)**3 + Ol)**0.5)

d_dm = -(0.5*c*(x + 1)*(z + 1)**3)/(H0*(Om*(z + 1)**3 + y)**1.5)

d_dy = 








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
    h  = Parameter(default=70.)
    Om = Parameter(default=0.3)
    Ol = Parameter(default=0.7)

    @staticmethod
    def evaluate(x, H0, Om, Ol):
        c = 2.99792458E+5
        return (c*(1.+x)/H0) * (((1.+x)**3.)*Om + Ol)**-0.5
    
    @staticmethod
    def fit_deriv(x, H0, Om, Ol):
        """
        h : H0
        m : Om
        L : Ol
        These are the derivatives of the DL func for the Lambda CDM model. 
        x - the z that is not part of the integral.
        """
        c = 2.99792458E+5
        d_H0  = -(c*(x + 1))/(h**2*(Om*(x + 1)**3 + Ol)**0.5)
        d_Om = -(0.5*c*(x + 1)**4)/(h*(Om*(x + 1)**3 + Ol)**1.5)
        d_Ol = -(0.5*c*(x + 1))/(h*(Om*(x + 1)**3 + Ol)**1.5)
        d_z  = 
        return [d_h, d_Om, d_Ol, d_z]



# FlatLambdaCDM 
redshifts = np.linspace(1, 8, 100)
cosmoconstants = {'H0':67.8, 'Om': 0.308}
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='FlatLambdaCDM', 
                                 units='Mpc',
                                 **cosmoconstants)

# x2 : redshifts,  y2 : LumDist
x2 = redshifts
y2 = DL
np.random.seed(1) #seed
y2 += np.random.normal(0., 0.1, x2.shape) # 0.5
y2_err = np.ones(x2.shape)*0.1


plt.clf()
plt.errorbar(x2, y2, yerr=y2_err, fmt='.')
plt.show()


m = FLCDM(h=71., Om=0.308, Ol=1-0.308)
fitter = fitting.LevMarLSQFitter()
f = fitter(model=m, x=x2, y=y2, weights=1.0/y2_err**2)

plt.errorbar(x2, y2, yerr=y2_err, fmt='.')
plt.plot(x2, f(x2))
plt.yscale('log')
plt.show()



plt.errorbar(x2, y2, yerr=y2_err, fmt='.')
plt.plot(x2,sine_fit(x2))
plt.show()

print(sine_fit)

calc_reduced_chi_square(sine_fit(x2), x2, y2, y2_err, len(x2), 3)


