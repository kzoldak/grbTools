from numpy import sqrt
from scipy import integrate


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
    Default cosmology is taken from the Planck 2015 results at
        https://arxiv.org/pdf/1502.01589.pdf
    These are:
    H_knot (hubble constant): 67.8 +- 

    

    """
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
            H_knot = [j for j in cc.keys() if j in hc_options][0]
            omega_M = [j for j in cc.keys() if j in md_options][0]
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
        eqn = (1./(sqrt(((1.+z)*(1.+z)*(1.+omega_m*z))-(z*(2.+z)*omega_l ))))
        return eqn
    # Integrate from 0 to z, [0] ignores the error on integration.
    DL = integrate.quad(_lumdist_equation, 0.0, z)[0]
    # Lum Dist in Mpc
    DL_Mpc  = (c*(1.+z)/H_knot) * DL  
    # Lum Dist in cm is what we want.
    DL_cm   = DL_Mpc * Mpctocm 
    return DL_cm



# # **********************************************************
         
#     if cc is not None and isinstance(cc, dict):
#         # Check for correct keywords in dict
#         hc_options = ['hubble_constant', 'H_knot']
#         if any(i in cc.keys() for i in hc_options):
#             pass
#         else:
#             msg = ("Must use '{}' or '{}' for "
#                    "Hubble Constant.".format(*hc_options))
#             raise Exception(msg)
#         # Check for correct keywords in dict
#         md_options = ['matter_density', 'omega_M']
#         if any(i in cc.keys() for i in md_options):
#             pass
#         else:
#             msg = ("Must use '{}' or '{}' for "
#                    "Matter Density.".format(*md_options))
#             raise Exception(msg)
        
#         H_knot = [j for j in cc.keys() if j in hc_options][0]
#         omega_M = [j for j in cc.keys() if j in md_options][0]



        
