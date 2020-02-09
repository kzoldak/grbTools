#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 10:22:52 2020

@author: kimzoldak
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 21:29:28 2020

@author: kimzoldak
"""




import numpy as np

from numpy import exp

import matplotlib.pyplot as plt
%matplotlib inline

import astropy
from astropy.modeling import Fittable1DModel, Parameter,m fitting


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



class FlatLambdaCDM(FIttable1dModel):
    """
    All parameters go here.
    """
    #c = Parameter(default=1.)
    #H0 = Parameter(default=1.)
    Om = Parameter(default=0.27)
    Ol = Parameter(default=1.-0.27)
    z = Parameter(default=1.)

    @staticmethod
    def evaluate(x, z, Om, Ol):
        H0 = 71.
        c = 2.99792458E+5
        return (c*(1.+x)/H0) * (((1.+z)**3.)*Om + Ol)**-0.5
    
    @staticmethod
    def fit_deriv(x, z, Om, Ol):
        """
        These are the derivatives of the DL func for the Lambda CDM model. 
        x - the z that is not part of the integral.
        """
        H0 = 71.
        c = 2.99792458E+5
        d_z = (-1.5*Om*(1+z)**2)/((Om*(1+z)**3 + Ol)**1.5)
        d_Om = -(0.5*(1+z)**3)/(Om*(1+z)**3 + Ol)**1.5
        d_Ol = -0.5/(Om*(1+z)**3 + Ol)**1.5
        A = (c*(1.+x)/H0)
        d_z = A * d_z
        d_Om = A * d_Om
        d_Ol = A * d_Ol
        # np.ones(x.shape)
        return [d_z, d_Om, d_Ol]
    



#x2 = np.linspace(0,10,100)
#a = 3
#b = 2
#c = 4
#d = 1
#y2 = a*np.sin(b*x2+c)+d
#np.random.seed(1) #seed
#y2 += np.random.normal(0., 0.5, x2.shape)
#y2_err = np.ones(x2.shape)*0.3


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
#plt.yscale('log')
plt.show()


# (x, z, Om, Ol)
model = FlatLambdaCDM(a=4.,b=2.,c=4.,d=0.)
fitter = fitting.LevMarLSQFitter()
sine_fit = fitter(sine_model, x2, y2, weights = 1.0/y2_err**2)






sine_model = SineNew(a=4.,b=2.,c=4.,d=0.)
fitter = fitting.LevMarLSQFitter()
sine_fit = fitter(sine_model, x2, y2, weights = 1.0/y2_err**2)



plt.errorbar(x2, y2, yerr=y2_err, fmt='.')
plt.plot(x2,sine_fit(x2))
plt.show()

print(sine_fit)

calc_reduced_chi_square(sine_fit(x2), x2, y2, y2_err, len(x2), 3)


