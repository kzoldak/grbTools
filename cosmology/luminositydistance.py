"""
At the end of this file, we provide a very basic LumDist function that 
I used so much that I can't get used to not using it. It assumes the 
Flat Lambda CDM cosmology and it's all set up they way I need it for 
calculating GRB luminosity distances. 

"""

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
    # def __init__(self, **constants):
    #     self.models = ['FlatLambdaCDM','w0wpCDM','w0wpCDM_error',
    #                    'w0waCDM', 'weylgravity']
    #     if constants:
    #         self.H0 = constants['H0']
    #         try:
    #             self.Om = constants['Om']
    #             self.Ol = 1.0 - self.Om
    #         except:
    #             pass
    #     else:
    #         self.H0 = 67.8
    #         self.Om = 0.308
    #         self.Ol = 1.0 - self.Om

    def __init__(self, H0=None, Om=None, Ol=None, w0=None, wp=None, 
                    wa=None, q0=None):
        self.models = ['FlatLambdaCDM','w0wpCDM','w0wpCDM_error',
                       'w0waCDM', 'weylgravity']
        if H0 is not None:
            self.H0 = H0
        if Om is not None:
            self.Om = Om
        if Ol is not None:
            self.Ol = Ol
        if w0 is not None:
            self.w0 = w0
        if wp is not None:
            self.wp = wp
        if wa is not None:
            self.wa = wa
        if q0 is not None:
            self.q0 = q0

    # *** NOTES ON UNIT CONVERSIONS ***
    # Make `_pc_to_cm` your home base. Then take that value and convert it to 
    #  other SI units from there using e.g., _cm_to_km or _cm_to_m, etc.
    #  We do this because these two units are used in our luminosity distance
    #  calculations. DLs in cm are used for computing Eiso energies and DLs in 
    #  parsecs is used for computing mu magnitudes in the distance modulus 
    #  equation. 
    # Units below need to be finished later. 
    
    @property
    def _Mpc_to_pc(self):
        # 1 Mpc = 1E6 pc
        return 1.0E6
    @property
    def _pc_to_Mpc(self):
        # 1 Mpc = 1E6 pc
        return 1.0E-6
    @property
    def _pc_to_cm(self):
        # 1 pc = 3.086e+18 cm
        return 3.08567758128E+18
    @property
    def _km_to_cm(self):
        # 1 km = 100000 cm
        return 1.0E5
    @property
    def _cm_to_km(self):
        # 1 km = 100000 cm
        return 1.0E-5
    
    @property
    def _speedoflight(self):
        """
        SPEED OF LIGHT in  Units: km/s
        If you want other units, convert it later. 
        E.g.,
        c = self._speedoflight * self._km_to_cm
        if you want c in cm/s instead of km/s
        """
        c = 2.99792458E+5 # km/s
        return c

    def _lumdist_FlatLambdaCDM(self, redshift):
        """
        This is the function we use in our work, but different cosmo constants
        as well as DL units.
        """
        z = redshift
        H0 = self.H0
        Om = self.Om
        try:
            Ol = self.Ol
        except:
            Ol = 1.0 - Om
        c = self._speedoflight
        pt0 = (c*(1.+z)/H0)
        def Aint(z):
            pt1 = ((1.+z)**3.)*Om
            pt2 = Ol
            return (pt1 + pt2)**-0.5
        AA = integrate.quad(Aint, 0.0, z)
        dl = pt0 * AA[0]
        return dl


    def _lumdist_w0wpCDM(self, redshift):
        """
        Lower order expansion for dark energy term (w). 
        Equation 14 in Riess et al. 2004.
        w(z) = w0 + w'z
        This is what I used to call Riess Cosmo. 
        """
        z = redshift
        H0 = self.H0
        Om = self.Om
        try:
            Ol = self.Ol
        except:
            Ol = 1.0 - Om
        try: 
            w0 = self.w0
        except:
            w0 = -1.31 
        try:
            wp = self.wp
        except:
            wp = 1.48 
        c = self._speedoflight
        pt0 = (c*(1.+z)/H0)
        def Aint(z):
            pt1 = ((1.+z)**3.)*Om
            pt2 = Ol
            pt3 = ((1.+z)**(3*(1.+w0-wp)))*np.exp(3.*wp*z)
            return (pt1 + pt2 * pt3)**-0.5
        AA = integrate.quad(Aint, 0.0, z)
        dl = pt0 * AA[0]
        return dl 


    def _lumdist_w0wpCDM_error(self, redshift):
        """
        Lower order expansion for dark energy term (w). 
        Equation 14 in Riess et al. 2004.
        w(z) = w0 + w'z
        """
        z = redshift
        H0 = self.H0
        Om = self.Om
        try:
            Ol = self.Ol
        except:
            Ol = 1.0 - Om
        try: 
            w0 = self.w0
        except:
            w0 = -1.31 
        try:
            wp = self.wp
        except:
            wp = 1.48 
        c = self._speedoflight
        pt0 = (c*(1.+z)/H0)
        def Aint(z):
            pt1 = ((1.+z)**3.)*Om
            pt2 = Ol
            pt3 = ((1.+z)**(3*(1.+w0-wp)))*np.exp(-3.*wp*z)
            return (pt1 + pt2 * pt3)**-0.5
        AA = integrate.quad(Aint, 0.0, z)
        dl = pt0 * AA[0]
        return dl 


    def _lumdist_w0waCDM(self, redshift):
        """
        Linder 2003 the more complex expansion of w(z).
        w(a) = w_0 + w_a * (1.-a) = w_0 + w_a*z/(1.+z)
        or
        w(z) = w0 + wa*(z/(1.+z))
        DL eqn also given as 
        Lower order expansion for dark energy term (w). 
        Equation 14 in Riess et al. 2004.
        w(z) = w0 + w'z
        """
        z = redshift
        H0 = self.H0
        Om = self.Om
        try:
            Ol = self.Ol
        except:
            Ol = 1.0 - Om
        try: 
            w0 = self.w0
        except:
            w0 = -0.9
        try:
            wa = self.wa
        except:
            wa = 0.2
        c = self._speedoflight
        pt0 = (c*(1.+z)/H0)
        def Aint(z):
            pt1 = ((1.+z)**3.)*Om
            pt2 = Ol
            pt3 = ((1.+z)**(3*(1.+w0+wa)))*np.exp(-3.*wa*(z/(1.+z)))
            return (pt1 + pt2 * pt3)**-0.5
        AA = integrate.quad(Aint, 0.0, z)
        dl = pt0 * AA[0]
        return dl 


    def _lumdist_weylgravity(self, redshift):
        """
        Weyl Gravity.
        Equation 237 in Mannheim 2006 paper
        q_knot      = -0.37 or -0.2
        """
        z = redshift
        H0 = self.H0
        try:
            q0 = self.q0
        except:
            q0 = -0.37
        c = self._speedoflight
        pt0 = -c*((1.+z)**2)/(H0*q0)
        pt1 = (1. + q0 - (q0/((1+z)**2)))**0.5
        dl = pt0 * (1. - pt1)
        return dl


    def luminosity_distance(self, redshifts, model='FlatLambdaCDM', 
                            units=None, **constants):
        """
        luminosity_distance(redshifts, model='FlatLambdaCDM', units=None, 
                            **constants)
        
        Parameters: 
        -----------
        redshifts : float, list, or array. Of redshifts. 
        model : str. One of the following model names:
                'FlatLambdaCDM','w0wpCDM','w0wpCDM_error',
                'w0waCDM', 'weylgravity'
                ** Default is 'FlatLambdaCDM'
        units : str, units to return luminosity distances in. 
                Options are: 'pc', 'Mpc', 'cm', None
                Default is None, which will return Mpc by default.
                If using for Eiso calculations, use units='cm'.
                If using for distance modulus, use units='pc' or 'Mpc'. 
        constants : dict, of cosmology constants associated with a model. 
              'FlatLambdaCDM'  : 'H0', 'Om', 'Ol'
              'w0wpCDM'        : 'H0', 'Om', 'Ol', 'w0', 'wp'
                                  w0 = -1.31, wp = 1.48
              'w0waCDM'        : 'H0', 'Om', 'Ol', 'w0', 'wa'
                                  w0 = -0.9, wa = 0.2
              'weylgravity'    : 'H0', 'q0'
                                  q0 = -0.37
                There are no defaults for H0 or Om for any model, these must 
                be specified. For all cases here, Ol = 1-Om.
                For the rest of the constants, if a value isn't specified, the 
                default will be used; defaults listed below each model above. 
        
        Notes: 
        ------
        'FlatLambdaCDM' is the concordance model (most accepted). 
          Hogg 2000, 
            “Distance measures in cosmology”, 
            https://arxiv.org/abs/astro-ph/9905116
          Schaefer 2007 (eqn 10), 
            "The Hubble Diagram to Redshift >6 from 69 Gamma-Ray Bursts"
            https://iopscience.iop.org/article/10.1086/511742   or   
            https://arxiv.org/abs/astro-ph/0612285
        
        'w0wpCDM' is the w(z) expansion.
          Linder 2003 (pg 1, Sec A. Linear w(z)), 
            “Exploring the Expansion History of the Universe”, 
            https://arxiv.org/abs/astro-ph/0208512
          Riess et al. 2004 (eqn 14), 
            "Type Ia Supernova Discoveries at z>1 From the Hubble Space 
            Telescope: Evidence for Past Deceleration and 
            Constraints on Dark Energy Evolution"
            https://iopscience.iop.org/article/10.1086/383612   or   
            https://arxiv.org/abs/astro-ph/0402512
            ** w0 and wp values adopted from Riess and Schaefer. Riess used SN.
        
        'w0waCDM' is the w(z) = w0 + wa * z/(1+z) expansion.
          Linder 2003 (pg 2, B. A new parametrization of the equation of state),
            “Exploring the Expansion History of the Universe”, 
            https://arxiv.org/abs/astro-ph/0208512
          Schaefer 2007 (eqn 11), 
            "The Hubble Diagram to Redshift >6 from 69 Gamma-Ray Bursts"
            https://iopscience.iop.org/article/10.1086/511742   or   
            https://arxiv.org/abs/astro-ph/0612285
            ** w0 and wa values adopted from Schaefer, who used GRBs and their 
            luminosity distances relations.
            
        'weylgravity' 
        Mannheim 2006 (eqn 237).
        https://arxiv.org/abs/astro-ph/0505266
        
        
        """
        redshifts = np.asarray(redshifts)
        if constants:
            for key in constants.keys():
                self.__setattr__(key, constants[key])

        self.models = ['FlatLambdaCDM','w0wpCDM','w0wpCDM_error',
                       'w0waCDM', 'weylgravity']
        if model not in self.models:
            raise Exception('Model not an option. Options are: %r'%self.models)
        
        if model == 'FlatLambdaCDM':
            func = self._lumdist_FlatLambdaCDM
        elif model == 'w0wpCDM':
            func = self._lumdist_w0wpCDM
        elif model == 'w0wpCDM_error':
            func = self._lumdist_w0wpCDM_error
        elif model == 'w0waCDM':
            func = self._lumdist_w0waCDM
        elif model == 'weylgravity':
            func = self._lumdist_weylgravity
        else:
            pass
        
        # Send of Calculation for Luminosity Distance.
        distances = np.asarray([func(z) for z in redshifts])

        if units is not None:
            self.units = units
        unit_lst = ['pc', 'Mpc', 'cm']
        if (units is not None) and (units not in unit_lst): 
            #raise Exception("Currently only recognize one of following: %r "%unit_list)
            raise Exception("Currently don't recognize units.")
        if units == 'cm':
            # Mpc to pc  then  pc to cm
            return distances * self._Mpc_to_pc * self._pc_to_cm
        elif units == 'pc':
            return distances * self._Mpc_to_pc  # convert Mpc to pc
        else:
            return distances # return default in Mpc
        



        # Derivatives of the lumionosity distance functions above.
        # THIS IS NOT FINISHED!!! 
        def _deriv_wrt_z__FlatLCDM_lumdist(z, Om, Ol):
            return (-1.5*Om*(1+z)**2)/((Om*(1+z)**3 + Ol)**1.5)
        
        
        def _deriv_wrt_z__w0wp_lumdist(z, Om, Ol, w0, wp):
            from numpy import exp
            x = Om
            y = Ol
            a = wp
            b = w0
            return -(0.5*(3*a*y*exp(3*a*z)*(1+z)**(3*(-a+b+1))+3*y*(-a+b+1) \
                          *exp(3*a*z)*(1+z)**(3*(-a+b+1)-1)+3*x*(1+z)**2))/ \
                          (y*e**(3*a*z)*(1+z)**(3*(-a+b+1))+x*(1+z)**3)**1.5
        
        
        def _deriv_wrt_z__w0wa_lumdist(z, Om, Ol, w0, wa):
            from numpy import exp
            y = Ol
            x = Om
            a = wa
            b = w0
            return -(0.5*(y*exp(-(3*a*z)/(1+z))*((3*a*z)/(1+z)**2 -(3*a)/ \
                        (1+z))*(1+z)**(3*(a+b+1))+3*y*(a+b+1)*exp(-(3*a*z)/ \
                        (1+z))*(1+z)**(3(a+b+1)-1)+3*x*(1+z)**2))/(y* \
                        exp(-(3*a*z)/(1+z))*(1+z)**(3*(a+b+1))+x*(1+z)**3)**1.5






def LumDist(redshift, cosmoconstants=None):
    """
    Parameters:
    -----------
    redshift: float, redshift of the object.
    cosmoconstants: list or dict of cosmology constants.
                    Only provide the hubble constant and the 
                    matter density. Energy density is computed from
                    the matter density:   
                    omega_L = 1 - omega_M
                    where L stands from lambda and M stands for matter.

                    Default values are:
                    H_knot = 67.8
                    omega_M = 0.308
                    omega_L = 1 - omega_M
                    ** See Notes for more details! **

                    If you want to pass your own constants, 
                    you can pass them as a list or a dictionary. 
                    If using a list, order matters. If using a 
                    dictionary, order is irrelevant.

                    List structure:
                          [hubble_constant, matter_density]
                            or as better known as:
                          [H_knot, omega_M]
                    
                    Dict structure:
                        {'hubble_constant': value,
                         'matter_density': value}
                              or 
                         {'H_knot': value,
                          'omega_M': value}

                    All other parameter names will throw and error. 

    Notes:
    ------

    DEFAULT COSMOLOGY:
    The default Cosmology Constants are:
    H_knot      = 67.8     +- 0.09  km s^-1 Mpc^-1    Hubble Constant
    omega_M     = 0.308    +- 0.012                   Matter Density
    omega_L     = 1.0 - omega_M                       Energy Density
    These are from the Planck 2015 results at 
        https://arxiv.org/pdf/1502.01589.pdf

    OUTDATED COSMOLOGY:
    Outdated Cosmology Parameters that were once used in my analyses.
    MAY HAVE USED THESE FOR RMFIT.
    H_knot = 67.3
    omega_M = 0.315
    omega_L = 1 - omega_M
    Was originally used on GRBs modeled with RMFIT.
    # ***************
    H_knot      = 70.0
    omega_M     = 0.27
    omega_L     = 1.0 - omega_M
    XSPEC and BXA results orig used these constants. 

    ** Important Note: The Eiso energies computed in all of my 
        GRB analyses are now consistent in cosmology and use the 
        default parameters. **

    
    We chose this cosmology because 
    2nd GBM Spectral Catalog states:  
    ' All current Cosmological parameters are obtained from the 
      Planck Collaboration: Ade et al. 2013 '

    """
    #print(cosmoconstants)
    cc = cosmoconstants
    if cc is not None:
        if isinstance(cc, list):  # if list format
            # Assign 
            H_knot, omega_M = cc
        # It is better to use isinstance(cc, dict) to check type instead
        # of using if type(cc) is dict or type(cc) == dict because
        # neither of those will catch the difference between traditional
        # dictionaries and OrderedDict dictionaries.
        #   If not sure, see documentation for it after importing.
        #   from collections import OrderedDict
        elif isinstance(cc, dict):  # if dict format
            # Check for correct keywords in dict
            hc_options = ['hubble_constant', 'H_knot']
            if any(i in cc.keys() for i in hc_options):
                pass
            else:
                msg = ("Must use '{}' or '{}' for "
                       "Hubble Constant.".format(*hc_options))
                raise Exception(msg)
            # Check for correct keywords in dict
            md_options = ['matter_density', 'omega_M']
            if any(i in cc.keys() for i in md_options):
                pass
            else:
                msg = ("Must use '{}' or '{}' for "
                       "Matter Density.".format(*md_options))
                raise Exception(msg)
            # Assign constants
            # Everything inside cc[] is the keyword we want for cc
            # i.e., below is same as cc['H_knot'] 
            #   or cc['hubble_constant'], whichever key word was used.
            H_knot = cc[ [j for j in cc.keys() if j in hc_options][0] ]
            omega_M = cc[ [j for j in cc.keys() if j in md_options][0] ]
        else:
            dtype = type(cosmoconstants)
            msg = ("cosmoconstants must be list or dict type, "
                   "got {dtype!r}").format(dtype=dtype)
            raise TypeError(msg)
    else: 
        # If no cosmo constants provided, use these as defaults:
        H_knot = 67.8
        omega_M = 0.308
    omega_L = 1.0 - omega_M

    # Other constants:
    c           = 2.99792458e5    # SPEED OF LIGHT (km/s)
    Mpctocm     = 3.08567758e24   # Conversion of Mpc to cm 
    keVtoerg    = 1.60217657e-9   # Conversion of keV to erg
    MeVtoerg    = 1.60217657e-6   # Conversoinof MeV to erg

    z = redshift  # GRB redshift

    # Luminosity Distance Calculation:
    def _lumdist_equation(z):
        eqn = (1./(sqrt(((1.+z)*(1.+z)*(1.+omega_M*z))-(z*(2.+z)*omega_L))))
        return eqn
    # Integrate from 0 to z, [0] ignores the error on integration.
    DL = integrate.quad(_lumdist_equation, 0.0, z)[0]
    # Lum Dist in Mpc
    DL_Mpc  = (c*(1.+z)/H_knot) * DL  
    # Lum Dist in cm is what we want.
    DL_cm   = DL_Mpc * Mpctocm 
    return DL_cm


