from __future__ import division

import numpy as np
from numpy import pi, sqrt

from scipy import integrate


class Cosmology(object):
    """
    examples:
    
    kim = Cosmology()
    kim.luminosity_distance(redshifts=redshifts, model='c', **{'H0':65., 'Om':0.3})
    
    kim = Cosmology(**{'H0':72., 'Om':0.3})
    kim.luminosity_distance(redshifts=redshifts, model='c')
    
    # Default cosmology will be used. 
    kim = Cosmology()
    kim.luminosity_distance(redshifts=redshifts, model='c')
    
    
    """
    def __init__(self, **constants):
        if constants:
            self.H0 = constants['H0']
            try:
                self.Om = constants['Om']
                self.Ol = 1.0 - self.Om
            except:
                pass
        else:
            self.H0 = 67.8
            self.Om = 0.308
            self.Ol = 1.0 - self.Om

    def _lumdist_lcdm(self, redshift):
        """
        This is the function we use in our work, but different cosmo constants
        as well as DL units.
        """
        z = redshift
        H0 = self.H0
        Om = self.Om
        Ol = 1.0 - Om
        c = 2.99792458e5      # SPEED OF LIGHT    Units: km/s
        def Aint(z):
            zp1 = 1.+z 
            return ((zp1**3.)*Om + Ol)**-0.5
            #return (1./(np.sqrt(((1.+z)*(1.+z)*(1.+Om*z))-(z*(2.+z)*Ol))))
        AA = integrate.quad(Aint, 0.0, z)
        dl = (c*(1.+z)/H0)*AA[0]  # DL here has units:  Mpc
        dl = dl * (1.e6) # convert Mpc to pc
        return dl # in parsecs


    def _lumdist_riess(self, redshift):
        """
        Lower order expansion for dark energy term (w). 
        Equation 14 in Riess et al. 2004.
        w(z) = w0 + w'z
        """
        z = redshift
        H0 = self.H0
        Om = self.Om
        Ol = 1.0 - Om
        c = 2.99792458e5      # SPEED OF LIGHT    Units: km/s
        w0 = -1.31
        wp = 1.48    # p: prime symbol
        def Aint(z):
            zp1 = 1. + z
            Iscl = (zp1**(3.*(1.+w0-wp))) * np.exp(3.*wp*z)
            return ((zp1**3.)*Om+Ol*Iscl)**-0.5
        # def Aint(z):
        #     return 1./np.sqrt(((1.+z)**3)*Om+ \
        #                      Ol*((1.+z)**(3*(1+w0-wp)))*np.exp(3*wp*z))
        AA = integrate.quad(Aint, 0.0, z)
        dl = (c*(1.+z)/H0)*AA[0]  # DL here has units:  Mpc
        dl = dl * (1.e6) # convert Mpc to pc
        return dl # in parsecs

    def _lumdist_riess_error(self, redshift):
        """
        Lower order expansion for dark energy term (w). 
        Equation 14 in Riess et al. 2004.
        w(z) = w0 + w'z
        """
        z = redshift
        H0 = self.H0
        Om = self.Om
        Ol = 1.0 - Om
        c = 2.99792458e5      # SPEED OF LIGHT    Units: km/s
        w0 = -1.31
        wp = 1.48    # p: prime symbol
        def Aint(z):
            zp1 = 1. + z
            Iscl = (zp1**(3.*(1.+w0-wp))) * np.exp(-3.*wp*z)
            return ((zp1**3.)*Om+Ol*Iscl)**-0.5
        AA = integrate.quad(Aint, 0.0, z)
        dl = (c*(1.+z)/H0)*AA[0]  # DL here has units:  Mpc
        dl = dl * (1.e6) # convert Mpc to pc
        return dl # in parsecs


    # def _lumdist_linder(self, redshift):
    #     """
    #     Linder 2003 expansion.
    #     w(a) = w_0 + w_a * (1.-a) = w_0 + w_a*z/(1.+z)
    #     or
    #     w(z) = w0 + wa*z*(1.+z)^-1
    #     Lower order expansion for dark energy term (w). 
    #     Equation 14 in Riess et al. 2004.
    #     w(z) = w0 + w'z
    #     """
    #     z = redshift
    #     H0 = self.H0
    #     Om = self.Om
    #     Ol = 1.0 - Om
    #     c = 2.99792458e5      # SPEED OF LIGHT    Units: km/s
    #     w0 = 
    #     wa = 
    #     def Aint(z):
    #         zp1 = 1. + z
    #         Iscl = (zp1**(3.*(1.+w0-wp))) * np.exp(3.*wp*z)
    #         return ((zp1**3.)*Om+Ol*Iscl)**-0.5
    #     # def Aint(z):
    #     #     return 1./np.sqrt(((1.+z)**3)*Om+ \
    #     #                      Ol*((1.+z)**(3*(1+w0-wp)))*np.exp(3*wp*z))
    #     AA = integrate.quad(Aint, 0.0, z)
    #     dl = (c*(1.+z)/H0)*AA[0]  # DL here has units:  Mpc
    #     dl = dl * (1.e6) # convert Mpc to pc
    #     return dl # in parsecs


    def _lumdist_weylgravity(self, redshift):
        """
        Weyl Gravity.
        Equation 237 in Mannheim 2006 paper
        q_knot      = -0.37 or -0.2
        """
        z = redshift
        H0 = self.H0
        c = 2.99792458e5      # SPEED OF LIGHT    Units: km/s
        q0 = -0.37
        dl  = (-c*((1.+z)**2)/(H0*q0)) * \
                    (1 - np.sqrt(1+q0-(q0/((1+z)**2))))
        dl = dl * (1.e6) # convert Mpc to pc
        return dl # in parsecs
    
    
    def luminosity_distance(self, redshifts, model, **constants):
        redshifts = np.asarray(redshifts)
        if constants:
            self.H0 = constants['H0']
            try:
                self.Om = constants['Om']
                self.Ol = 1.0-self.Om
            except:
                pass
        else:
            self.H0 = self.H0
            self.Om = self.Om
            self.Ol = 1.0-self.Om
        if (model.startswith('c')) or (model.startswith('l')):
            # c - concordance, l - Lambda CDM
            func = self._lumdist_lcdm 
        elif model.startswith('r'):
            if 'err' in model:
                func = self._lumdist_riess_error
            else:
                func = self._lumdist_riess 
        elif model.startswith('w'):  
            func = self._lumdist_weylgravity 
        
        distances = np.asarray([func(z) for z in redshifts])
        return distances






