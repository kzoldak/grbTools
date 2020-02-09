from __future__ import division
from collections import OrderedDict


def powerlaw(dpar=None):
    pars = OrderedDict([
        ('PhoIndex__1', 
            '-energy*energy**(-PhoIndex__1)*norm__2*log(energy)'),
        ('norm__2', 
            'energy*energy**(-PhoIndex__1)')
        ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def lpow(dpar=None):
    pars = OrderedDict([
        ('plIndex__1', 
            'energy*norm__2*(0.01*energy)**plIndex__1*log(0.01*energy)'),
        ('norm__2', 
            'energy*(0.01*energy)**plIndex__1')
        ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def bbody(dpar=None):
    """
    Notes:
    ------
    When integrating the blackbody function's partial derivative wrt kT,
    you will run into the error: 
        "math domain error" 
    when kT is small. For a GRB with z=0 and kT < 28.2 keV, we would 
    get this error if integrating up to 10 MeV. 
    Thus, we need to use a maximum limit. 

    In order to avoid the math domain error, you should limit the 
    integration of the bbody fn's partial derivative wrt kT up to 
    and energy of emax = 354.8 * kT. For every 1 keV, the maximum 
    allowed energy increases by 354.89 keV. 
    For example:
    kT     max energy
    --     ----------
    1       354.89
    2       709.78
    3       1064.67
    10      3548.89
    The maximum allowed energy is 354.89 * kT, but to be conservative, 
    best to use 354.8 (or even lower). 

    The bbody will no longer be contributing to the total energy flux 
    long before you hit the maximum energy. At this point it is only 
    the primary model that is increasing flux. 

    If this limit becomes a problem, try reducing the function you are 
    integrating to one where bbody is removed (i.e., grbm instead of 
        grbm+bbody or grbm+lpow instead of grbm+bbody+lpow). Since the 
        bbody is not contributing at those higher energies, this is 
        acceptable to do. 


    """
    # The kt__1 value looks funny because it's separated onto 3 lines. 
    #  In Python, you can separate a string like this without using \
    #  if you surround it by (). It can span multiple lines without 
    #  causing issues. It is read as all one line of code. 
    pars = OrderedDict([
        ('kT__1',
            ('8.0525*energy**4*norm__2*exp(energy/kT__1)/'
             '(kT__1**6*(exp(energy/kT__1) - 1)**2) - 32.21*energy**3*'
             'norm__2/(kT__1**5*(exp(energy/kT__1) - 1))')
            ),
        ('norm__2', 
            '8.0525*energy**3/(kT__1**4*(exp(energy/kT__1) - 1))')
        ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars

def cutoffpl(dpar=None):
    pars = OrderedDict([
        ('PhoIndex__1',
            ('-energy*energy**(-PhoIndex__1)*norm__3*'
             'exp(-energy/HighECut__2)*log(energy)')
            ),
        ('HighECut__2',
            ('energy**2*energy**(-PhoIndex__1)*norm__3*'
             'exp(-energy/HighECut__2)/HighECut__2**2')),
        ('norm__3', 
            'energy*energy**(-PhoIndex__1)*exp(-energy/HighECut__2)')
        ])
    
    if dpar is not None:
        return pars[dpar]
    else:
        return pars

def grbm(dpar=None):
    pars = OrderedDict([
        ('alpha__1',
            [('energy*norm__4*(0.01*energy)**alpha__1*'
                'exp(-energy/tem__3)*log(0.01*energy)'),
             ('energy*norm__4*(0.01*energy)**beta__2*(0.01*tem__3*'
                '(alpha__1 - beta__2))**(alpha__1 - beta__2)*'
                '(log(0.01*tem__3*(alpha__1 - beta__2)) + 1.0)*'
                'exp(-alpha__1 + beta__2) - energy*norm__4*'
                '(0.01*energy)**beta__2*(0.01*tem__3*'
                '(alpha__1 - beta__2))**(alpha__1 - beta__2)*'
                'exp(-alpha__1 + beta__2)')
             ]),
        ('beta__2',
            ['0',
             ('energy*norm__4*(0.01*energy)**beta__2*(0.01*tem__3*'
                '(alpha__1 - beta__2))**(alpha__1 - beta__2)*'
                '(-log(0.01*tem__3*(alpha__1 - beta__2)) - 1.0)*'
                'exp(-alpha__1 + beta__2) + energy*norm__4*'
                '(0.01*energy)**beta__2*(0.01*tem__3*'
                '(alpha__1 - beta__2))**(alpha__1 - beta__2)*'
                'exp(-alpha__1 + beta__2)*log(0.01*energy) + '
                'energy*norm__4*(0.01*energy)**beta__2*(0.01*tem__3*'
                '(alpha__1 - beta__2))**(alpha__1 - beta__2)*'
                'exp(-alpha__1 + beta__2)')
             ]),
        ('tem__3',
            [('energy**2*norm__4*(0.01*energy)**alpha__1*'
                'exp(-energy/tem__3)/tem__3**2'),
             ('100.0*energy*norm__4*(0.01*energy)**beta__2*'
                '(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*'
                '(0.01*alpha__1 - 0.01*beta__2)*'
                'exp(-alpha__1 + beta__2)/tem__3')]),
        ('norm__4',
            ['energy*(0.01*energy)**alpha__1*exp(-energy/tem__3)',
             ('energy*(0.01*energy)**beta__2*(0.01*tem__3*'
                '(alpha__1 - beta__2))**(alpha__1 - beta__2)*'
                'exp(-alpha__1 + beta__2)')
             ])])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def sbpl(dpar=None):
    pars = OrderedDict([
        ('alpha__1',
            ('10**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/'
                'ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/'
                'ebreak__3)**(3.33333333333333/log(10))) + (-0.15*'
                'alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**'
                '(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**'
                '(3.33333333333333/log(10))))*energy*norm__4*(0.01*energy)**'
                '(0.5*alpha__1 + 0.5*beta__2)*(0.15*log(0.5*(100.0/ebreak__3)'
                '**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**'
                '(3.33333333333333/log(10))) - 0.15*log(0.5*(energy/ebreak__3)'
                '**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)'
                '**(3.33333333333333/log(10))))*log(10) + 0.5*10**(-(-0.15*'
                'alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**'
                '(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**'
                '(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*'
                'beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/'
                'log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/'
                'log(10))))*energy*norm__4*(0.01*energy)**(0.5*alpha__1 '
                '+ 0.5*beta__2)*log(0.01*energy)')),
        ('beta__2',
            ('10**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)'
                '**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**'
                '(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)'
                '*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + '
                '0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy*'
                'norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(-0.15*'
                'log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5'
                '*(100.0/ebreak__3)**(3.33333333333333/log(10))) + 0.15*'
                'log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + '
                '0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*log(10) '
                '+ 0.5*10**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/'
                'ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/'
                'ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + '
                '0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/'
                'log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/'
                'log(10))))*energy*norm__4*(0.01*energy)**(0.5*alpha__1 + '
                '0.5*beta__2)*log(0.01*energy)')),
        ('ebreak__3',
            ('10**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)'
                '**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**'
                '(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)'
                '*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + '
                '0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy'
                '*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*((-0.15*'
                'alpha__1 + 0.15*beta__2)*(1.66666666666667*(energy/ebreak__3)'
                '**(-3.33333333333333/log(10))/(ebreak__3*log(10)) - '
                '1.66666666666667*(energy/ebreak__3)**(3.33333333333333/'
                'log(10))/(ebreak__3*log(10)))/(0.5*(energy/ebreak__3)**'
                '(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**'
                '(3.33333333333333/log(10))) - (-0.15*alpha__1 + '
                '0.15*beta__2)*(1.66666666666667*(100.0/ebreak__3)**'
                '(-3.33333333333333/log(10))/(ebreak__3*log(10)) - '
                '1.66666666666667*(100.0/ebreak__3)**(3.33333333333333/'
                'log(10))/(ebreak__3*log(10)))/(0.5*(100.0/ebreak__3)**'
                '(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**'
                '(3.33333333333333/log(10))))*log(10)')),
        ('norm__4',
            ('10**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/'
                'ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/'
                'ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + '
                '0.15*beta__2)*log(0.5*(energy/ebreak__3)**'
                '(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**'
                '(3.33333333333333/log(10))))*energy*(0.01*energy)**(0.5*'
                'alpha__1 + 0.5*beta__2)'))
        ])
    
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def grbm_lpow(dpar=None):
    pars = OrderedDict([
        ('alpha__1',
            [('energy*norm__4*(0.01*energy)**alpha__1*exp(-energy/tem__3)'
                '*log(0.01*energy)'),
             ('energy*(norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 '
                '- beta__2))**(alpha__1 - beta__2)*(log(0.01*tem__3*(alpha__1 '
                '- beta__2)) + 1.0)*exp(-alpha__1 + beta__2) - norm__4*'
                '(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**'
                '(alpha__1 - beta__2)*exp(-alpha__1 + beta__2))')
             ]),
        ('beta__2',
            ['0',
            ('energy*(norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 '
                '- beta__2))**(alpha__1 - beta__2)*(-log(0.01*tem__3*(alpha__1 '
                '- beta__2)) - 1.0)*exp(-alpha__1 + beta__2) + norm__4*'
                '(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**'
                '(alpha__1 - beta__2)*exp(-alpha__1 + beta__2)*log(0.01*energy) '
                '+ norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))'
                '**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2))')
            ]),
         ('tem__3',
            [('energy**2*norm__4*(0.01*energy)**alpha__1*'
                'exp(-energy/tem__3)/tem__3**2'),
            ('100.0*energy*norm__4*(0.01*energy)**beta__2*(0.01*tem__3*'
                '(alpha__1 - beta__2))**(alpha__1 - beta__2)*(0.01*alpha__1 '
                '- 0.01*beta__2)*exp(-alpha__1 + beta__2)/tem__3')
            ]),
         ('norm__4',
            ['energy*(0.01*energy)**alpha__1*exp(-energy/tem__3)',
             ('energy*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 '
                '- beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2)')
             ]),
         ('plIndex__5',
            ['energy*norm__6*(0.01*energy)**plIndex__5*log(0.01*energy)',
            'energy*norm__6*(0.01*energy)**plIndex__5*log(0.01*energy)']),
         ('norm__6',
            ['energy*(0.01*energy)**plIndex__5', 
             'energy*(0.01*energy)**plIndex__5'])])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars

def grbm_powerlaw(dpar=None):
    pars = OrderedDict([
        ('alpha__1',
            ['energy*norm__4*(0.01*energy)**alpha__1*exp(-energy/tem__3)*log(0.01*energy)',
            'energy*(norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(log(0.01*tem__3*(alpha__1 - beta__2)) + 1.0)*exp(-alpha__1 + beta__2) - norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2))']),
         ('beta__2',
            ['0',
            'energy*(norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(-log(0.01*tem__3*(alpha__1 - beta__2)) - 1.0)*exp(-alpha__1 + beta__2) + norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2)*log(0.01*energy) + norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2))']),
         ('tem__3',
            ['energy**2*norm__4*(0.01*energy)**alpha__1*exp(-energy/tem__3)/tem__3**2',
            '100.0*energy*norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(0.01*alpha__1 - 0.01*beta__2)*exp(-alpha__1 + beta__2)/tem__3']),
         ('norm__4',
            ['energy*(0.01*energy)**alpha__1*exp(-energy/tem__3)',
            'energy*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2)']),
         ('PhoIndex__5',
            ['-energy*energy**(-PhoIndex__5)*norm__6*log(energy)',
            '-energy*energy**(-PhoIndex__5)*norm__6*log(energy)']),
         ('norm__6',
            ['energy*energy**(-PhoIndex__5)', 
            'energy*energy**(-PhoIndex__5)'])
         ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def grbm_bbody(dpar=None):
    '''
    see bbody.__doc__
    '''
    pars = OrderedDict([
        ('alpha__1',
            ['energy*norm__4*(0.01*energy)**alpha__1*exp(-energy/tem__3)*log(0.01*energy)',
             'energy*(norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(log(0.01*tem__3*(alpha__1 - beta__2)) + 1.0)*exp(-alpha__1 + beta__2) - norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2))']),
         ('beta__2',
            ['0',
             'energy*(norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(-log(0.01*tem__3*(alpha__1 - beta__2)) - 1.0)*exp(-alpha__1 + beta__2) + norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2)*log(0.01*energy) + norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2))']),
         ('tem__3',
            ['energy**2*norm__4*(0.01*energy)**alpha__1*exp(-energy/tem__3)/tem__3**2',
             '100.0*energy*norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(0.01*alpha__1 - 0.01*beta__2)*exp(-alpha__1 + beta__2)/tem__3']),
         ('norm__4',
            ['energy*(0.01*energy)**alpha__1*exp(-energy/tem__3)',
             'energy*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2)']),
         ('kT__5',
            ['energy*(8.0525*energy**3*norm__6*exp(energy/kT__5)/(kT__5**6*(exp(energy/kT__5) - 1)**2) - 32.21*energy**2*norm__6/(kT__5**5*(exp(energy/kT__5) - 1)))',
             'energy*(8.0525*energy**3*norm__6*exp(energy/kT__5)/(kT__5**6*(exp(energy/kT__5) - 1)**2) - 32.21*energy**2*norm__6/(kT__5**5*(exp(energy/kT__5) - 1)))']),
         ('norm__6',
            ['8.0525*energy**3/(kT__5**4*(exp(energy/kT__5) - 1))',
             '8.0525*energy**3/(kT__5**4*(exp(energy/kT__5) - 1))'])
         ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars

def grbm_bbody_lpow(dpar=None):
    '''
    see bbody.__doc__
    '''
    pars = OrderedDict([
        ('alpha__1',
            ['energy*norm__4*(0.01*energy)**alpha__1*exp(-energy/tem__3)*log(0.01*energy)',
             'energy*(norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(log(0.01*tem__3*(alpha__1 - beta__2)) + 1.0)*exp(-alpha__1 + beta__2) - norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2))']),
         ('beta__2',
            ['0',
             'energy*(norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(-log(0.01*tem__3*(alpha__1 - beta__2)) - 1.0)*exp(-alpha__1 + beta__2) + norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2)*log(0.01*energy) + norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2))']),
         ('tem__3',
            ['energy**2*norm__4*(0.01*energy)**alpha__1*exp(-energy/tem__3)/tem__3**2',
             '100.0*energy*norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(0.01*alpha__1 - 0.01*beta__2)*exp(-alpha__1 + beta__2)/tem__3']),
         ('norm__4',
            ['energy*(0.01*energy)**alpha__1*exp(-energy/tem__3)',
             'energy*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2)']),
         ('kT__5',
            ['energy*(8.0525*energy**3*norm__6*exp(energy/kT__5)/(kT__5**6*(exp(energy/kT__5) - 1)**2) - 32.21*energy**2*norm__6/(kT__5**5*(exp(energy/kT__5) - 1)))',
             'energy*(8.0525*energy**3*norm__6*exp(energy/kT__5)/(kT__5**6*(exp(energy/kT__5) - 1)**2) - 32.21*energy**2*norm__6/(kT__5**5*(exp(energy/kT__5) - 1)))']),
         ('norm__6',
            ['8.0525*energy**3/(kT__5**4*(exp(energy/kT__5) - 1))',
             '8.0525*energy**3/(kT__5**4*(exp(energy/kT__5) - 1))']),
         ('plIndex__7',
            ['energy*norm__8*(0.01*energy)**plIndex__7*log(0.01*energy)',
             'energy*norm__8*(0.01*energy)**plIndex__7*log(0.01*energy)']),
         ('norm__8',
            ['energy*(0.01*energy)**plIndex__7', 'energy*(0.01*energy)**plIndex__7'])
         ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars

def grbm_bbody_powerlaw(dpar=None):
    '''
    see bbody.__doc__
    '''
    pars = OrderedDict([
        ('alpha__1',
            ['energy*norm__4*(0.01*energy)**alpha__1*exp(-energy/tem__3)*log(0.01*energy)',
            'energy*(norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(log(0.01*tem__3*(alpha__1 - beta__2)) + 1.0)*exp(-alpha__1 + beta__2) - norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2))']),
        ('beta__2',
          ['0',
           'energy*(norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(-log(0.01*tem__3*(alpha__1 - beta__2)) - 1.0)*exp(-alpha__1 + beta__2) + norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2)*log(0.01*energy) + norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2))']),
        ('tem__3',
          ['energy**2*norm__4*(0.01*energy)**alpha__1*exp(-energy/tem__3)/tem__3**2',
           '100.0*energy*norm__4*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*(0.01*alpha__1 - 0.01*beta__2)*exp(-alpha__1 + beta__2)/tem__3']),
        ('norm__4',
          ['energy*(0.01*energy)**alpha__1*exp(-energy/tem__3)',
           'energy*(0.01*energy)**beta__2*(0.01*tem__3*(alpha__1 - beta__2))**(alpha__1 - beta__2)*exp(-alpha__1 + beta__2)']),
        ('kT__5',
          ['energy*(8.0525*energy**3*norm__6*exp(energy/kT__5)/(kT__5**6*(exp(energy/kT__5) - 1)**2) - 32.21*energy**2*norm__6/(kT__5**5*(exp(energy/kT__5) - 1)))',
           'energy*(8.0525*energy**3*norm__6*exp(energy/kT__5)/(kT__5**6*(exp(energy/kT__5) - 1)**2) - 32.21*energy**2*norm__6/(kT__5**5*(exp(energy/kT__5) - 1)))']),
        ('norm__6',
          ['8.0525*energy**3/(kT__5**4*(exp(energy/kT__5) - 1))',
           '8.0525*energy**3/(kT__5**4*(exp(energy/kT__5) - 1))']),
        ('PhoIndex__7',
          ['-energy*energy**(-PhoIndex__7)*norm__8*log(energy)',
           '-energy*energy**(-PhoIndex__7)*norm__8*log(energy)']),
        ('norm__8',
          ['energy*energy**(-PhoIndex__7)', 'energy*energy**(-PhoIndex__7)'])
         ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def cutoffpl_lpow(dpar=None):
    pars = OrderedDict([
        ('PhoIndex__1',
            '-energy*energy**(-PhoIndex__1)*norm__3*exp(-energy/HighECut__2)*log(energy)'),
         ('HighECut__2',
            'energy**2*energy**(-PhoIndex__1)*norm__3*exp(-energy/HighECut__2)/HighECut__2**2'),
         ('norm__3', 
            'energy*energy**(-PhoIndex__1)*exp(-energy/HighECut__2)'),
         ('plIndex__4', 
            'energy*norm__5*(0.01*energy)**plIndex__4*log(0.01*energy)'),
         ('norm__5', 
            'energy*(0.01*energy)**plIndex__4')
         ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def cutoffpl_powerlaw(dpar=None):
    pars = OrderedDict([
        ('PhoIndex__1',
            '-energy*energy**(-PhoIndex__1)*norm__3*exp(-energy/HighECut__2)*log(energy)'),
        ('HighECut__2',
            'energy**2*energy**(-PhoIndex__1)*norm__3*exp(-energy/HighECut__2)/HighECut__2**2'),
        ('norm__3', 
            'energy*energy**(-PhoIndex__1)*exp(-energy/HighECut__2)'),
        ('PhoIndex__4', 
            '-energy*energy**(-PhoIndex__4)*norm__5*log(energy)'),
        ('norm__5', 
            'energy*energy**(-PhoIndex__4)')
        ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def cutoffpl_bbody(dpar=None):
    '''
    see bbody.__doc__
    '''
    pars = OrderedDict([
        ('PhoIndex__1',
            '-energy*energy**(-PhoIndex__1)*norm__3*exp(-energy/HighECut__2)*log(energy)'),
        ('HighECut__2',
            'energy**2*energy**(-PhoIndex__1)*norm__3*exp(-energy/HighECut__2)/HighECut__2**2'),
        ('norm__3', 
            'energy*energy**(-PhoIndex__1)*exp(-energy/HighECut__2)'),
        ('kT__4',
            'energy*(8.0525*energy**3*norm__5*exp(energy/kT__4)/(kT__4**6*(exp(energy/kT__4) - 1)**2) - 32.21*energy**2*norm__5/(kT__4**5*(exp(energy/kT__4) - 1)))'),
        ('norm__5', 
            '8.0525*energy**3/(kT__4**4*(exp(energy/kT__4) - 1))')
        ])
    
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def cutoffpl_bbody_lpow(dpar=None):
    '''
    see bbody.__doc__
    '''
    pars = OrderedDict([
        ('PhoIndex__1',
            '-energy*energy**(-PhoIndex__1)*norm__3*exp(-energy/HighECut__2)*log(energy)'),
        ('HighECut__2',
            'energy**2*energy**(-PhoIndex__1)*norm__3*exp(-energy/HighECut__2)/HighECut__2**2'),
        ('norm__3', 
            'energy*energy**(-PhoIndex__1)*exp(-energy/HighECut__2)'),
        ('kT__4',
            'energy*(8.0525*energy**3*norm__5*exp(energy/kT__4)/(kT__4**6*(exp(energy/kT__4) - 1)**2) - 32.21*energy**2*norm__5/(kT__4**5*(exp(energy/kT__4) - 1)))'),
        ('norm__5', 
            '8.0525*energy**3/(kT__4**4*(exp(energy/kT__4) - 1))'),
        ('plIndex__6', 
            'energy*norm__7*(0.01*energy)**plIndex__6*log(0.01*energy)'),
        ('norm__7', 
            'energy*(0.01*energy)**plIndex__6')
        ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def cutoffpl_bbody_powerlaw(dpar=None):
    '''
    see bbody.__doc__
    '''
    pars = OrderedDict([
        ('PhoIndex__1',
            '-energy*energy**(-PhoIndex__1)*norm__3*exp(-energy/HighECut__2)*log(energy)'),
        ('HighECut__2',
            'energy**2*energy**(-PhoIndex__1)*norm__3*exp(-energy/HighECut__2)/HighECut__2**2'),
        ('norm__3', 
            'energy*energy**(-PhoIndex__1)*exp(-energy/HighECut__2)'),
        ('kT__4',
            'energy*(8.0525*energy**3*norm__5*exp(energy/kT__4)/(kT__4**6*(exp(energy/kT__4) - 1)**2) - 32.21*energy**2*norm__5/(kT__4**5*(exp(energy/kT__4) - 1)))'),
        ('norm__5', 
            '8.0525*energy**3/(kT__4**4*(exp(energy/kT__4) - 1))'),
        ('PhoIndex__6', 
            '-energy*energy**(-PhoIndex__6)*norm__7*log(energy)'),
        ('norm__7', 
            'energy*energy**(-PhoIndex__6)')
        ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def sbpl_bbody(dpar=None):
    '''
    see bbody.__doc__
    '''
    pars = OrderedDict([
        ('alpha__1',
            'energy*(10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(0.345387763949107*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) - 0.345387763949107*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10)))) + 0.5*10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*log(0.01*energy))'),
        ('beta__2',
            'energy*(10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(-0.345387763949107*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + 0.345387763949107*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10)))) + 0.5*10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*log(0.01*energy))'),
        ('ebreak__3',
            '10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(2.30258509299405*(-0.15*alpha__1 + 0.15*beta__2)*(1.66666666666667*(energy/ebreak__3)**(-3.33333333333333/log(10))/(ebreak__3*log(10)) - 1.66666666666667*(energy/ebreak__3)**(3.33333333333333/log(10))/(ebreak__3*log(10)))/(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))) - 2.30258509299405*(-0.15*alpha__1 + 0.15*beta__2)*(1.66666666666667*(100.0/ebreak__3)**(-3.33333333333333/log(10))/(ebreak__3*log(10)) - 1.66666666666667*(100.0/ebreak__3)**(3.33333333333333/log(10))/(ebreak__3*log(10)))/(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))))'),
        ('norm__4',
            '10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)'),
        ('kT__5',
            'energy*(8.0525*energy**3*norm__6*exp(energy/kT__5)/(kT__5**6*(exp(energy/kT__5) - 1)**2) - 32.21*energy**2*norm__6/(kT__5**5*(exp(energy/kT__5) - 1)))'),
        ('norm__6', 
            '8.0525*energy**3/(kT__5**4*(exp(energy/kT__5) - 1))')
        ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def sbpl_lpow(dpar=None):
    pars = OrderedDict([
        ('alpha__1',
            'energy*(10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(0.345387763949107*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) - 0.345387763949107*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10)))) + 0.5*10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*log(0.01*energy))'),
        ('beta__2',
            'energy*(10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(-0.345387763949107*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + 0.345387763949107*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10)))) + 0.5*10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*log(0.01*energy))'),
        ('ebreak__3',
            '10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(2.30258509299405*(-0.15*alpha__1 + 0.15*beta__2)*(1.66666666666667*(energy/ebreak__3)**(-3.33333333333333/log(10))/(ebreak__3*log(10)) - 1.66666666666667*(energy/ebreak__3)**(3.33333333333333/log(10))/(ebreak__3*log(10)))/(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))) - 2.30258509299405*(-0.15*alpha__1 + 0.15*beta__2)*(1.66666666666667*(100.0/ebreak__3)**(-3.33333333333333/log(10))/(ebreak__3*log(10)) - 1.66666666666667*(100.0/ebreak__3)**(3.33333333333333/log(10))/(ebreak__3*log(10)))/(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))))'),
        ('norm__4',
            '10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)'),
        ('plIndex__5', 
            'energy*norm__6*(0.01*energy)**plIndex__5*log(0.01*energy)'),
        ('norm__6', 
            'energy*(0.01*energy)**plIndex__5')
        ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars


def sbpl_powerlaw(dpar=None):
    pars = OrderedDict([
        ('alpha__1',
            'energy*(10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(0.345387763949107*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) - 0.345387763949107*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10)))) + 0.5*10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*log(0.01*energy))'),
        ('beta__2',
            'energy*(10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(-0.345387763949107*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + 0.345387763949107*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10)))) + 0.5*10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*log(0.01*energy))'),
        ('ebreak__3',
            '10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(2.30258509299405*(-0.15*alpha__1 + 0.15*beta__2)*(1.66666666666667*(energy/ebreak__3)**(-3.33333333333333/log(10))/(ebreak__3*log(10)) - 1.66666666666667*(energy/ebreak__3)**(3.33333333333333/log(10))/(ebreak__3*log(10)))/(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))) - 2.30258509299405*(-0.15*alpha__1 + 0.15*beta__2)*(1.66666666666667*(100.0/ebreak__3)**(-3.33333333333333/log(10))/(ebreak__3*log(10)) - 1.66666666666667*(100.0/ebreak__3)**(3.33333333333333/log(10))/(ebreak__3*log(10)))/(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))))'),
        ('norm__4',
            '10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)'),
        ('PhoIndex__5', 
            '-energy*energy**(-PhoIndex__5)*norm__6*log(energy)'),
        ('norm__6', 
            'energy*energy**(-PhoIndex__5)')])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars



def sbpl_bbody_lpow(dpar=None):
    '''
    see bbody.__doc__
    '''
    pars = OrderedDict([
        ('alpha__1',
            'energy*(10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(0.345387763949107*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) - 0.345387763949107*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10)))) + 0.5*10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*log(0.01*energy))'),
        ('beta__2',
            'energy*(10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(-0.345387763949107*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + 0.345387763949107*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10)))) + 0.5*10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*log(0.01*energy))'),
        ('ebreak__3',
            '10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(2.30258509299405*(-0.15*alpha__1 + 0.15*beta__2)*(1.66666666666667*(energy/ebreak__3)**(-3.33333333333333/log(10))/(ebreak__3*log(10)) - 1.66666666666667*(energy/ebreak__3)**(3.33333333333333/log(10))/(ebreak__3*log(10)))/(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))) - 2.30258509299405*(-0.15*alpha__1 + 0.15*beta__2)*(1.66666666666667*(100.0/ebreak__3)**(-3.33333333333333/log(10))/(ebreak__3*log(10)) - 1.66666666666667*(100.0/ebreak__3)**(3.33333333333333/log(10))/(ebreak__3*log(10)))/(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))))'),
        ('norm__4',
            '10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)'),
        ('kT__5',
            'energy*(8.0525*energy**3*norm__6*exp(energy/kT__5)/(kT__5**6*(exp(energy/kT__5) - 1)**2) - 32.21*energy**2*norm__6/(kT__5**5*(exp(energy/kT__5) - 1)))'),
        ('norm__6', 
            '8.0525*energy**3/(kT__5**4*(exp(energy/kT__5) - 1))'),
        ('plIndex__7', 
            'energy*norm__8*(0.01*energy)**plIndex__7*log(0.01*energy)'),
        ('norm__8', 
            'energy*(0.01*energy)**plIndex__7')])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars



def sbpl_bbody_powerlaw(dpar=None):
    '''
    see bbody.__doc__
    '''
    pars = OrderedDict([
        ('alpha__1',
            'energy*(10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(0.345387763949107*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) - 0.345387763949107*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10)))) + 0.5*10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*log(0.01*energy))'),
        ('beta__2',
            'energy*(10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(-0.345387763949107*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + 0.345387763949107*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10)))) + 0.5*10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*log(0.01*energy))'),
        ('ebreak__3',
            '10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy*norm__4*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)*(2.30258509299405*(-0.15*alpha__1 + 0.15*beta__2)*(1.66666666666667*(energy/ebreak__3)**(-3.33333333333333/log(10))/(ebreak__3*log(10)) - 1.66666666666667*(energy/ebreak__3)**(3.33333333333333/log(10))/(ebreak__3*log(10)))/(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))) - 2.30258509299405*(-0.15*alpha__1 + 0.15*beta__2)*(1.66666666666667*(100.0/ebreak__3)**(-3.33333333333333/log(10))/(ebreak__3*log(10)) - 1.66666666666667*(100.0/ebreak__3)**(3.33333333333333/log(10))/(ebreak__3*log(10)))/(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))))'),
        ('norm__4',
            '10.0**(-(-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(100.0/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(100.0/ebreak__3)**(3.33333333333333/log(10))) + (-0.15*alpha__1 + 0.15*beta__2)*log(0.5*(energy/ebreak__3)**(-3.33333333333333/log(10)) + 0.5*(energy/ebreak__3)**(3.33333333333333/log(10))))*energy*(0.01*energy)**(0.5*alpha__1 + 0.5*beta__2)'),
        ('kT__5',
            'energy*(8.0525*energy**3*norm__6*exp(energy/kT__5)/(kT__5**6*(exp(energy/kT__5) - 1)**2) - 32.21*energy**2*norm__6/(kT__5**5*(exp(energy/kT__5) - 1)))'),
        ('norm__6', 
            '8.0525*energy**3/(kT__5**4*(exp(energy/kT__5) - 1))'),
        ('PhoIndex__7', 
            '-energy*energy**(-PhoIndex__7)*norm__8*log(energy)'),
        ('norm__8', 
            'energy*energy**(-PhoIndex__7)')
        ])
    if dpar is not None:
        return pars[dpar]
    else:
        return pars



def epeak_partials(model):
    """

    Parameters:
    -----------
    model: str, model name. 
            Options are: 'grbm', 'cutoffpl', 'sbpl'
            THESE ARE XSPEC ONLY!
    dpar: str or None, default is None.
            Pass None and get a dictionary of all the partial derivatives.
            Pass a single string and get only the partial derivative
              wrt that parameter. 

    Notes:
    ------
    The 'sbpl' model uses brkscale = 0.3. If you don't want to fix it 
    at this value, update this function and add brkscale to the parameters.
    We currently don't allow the function to use this parameter, so for 
    now it doesn't matter. 
    
    *********
    Propagating ebreak and tem errors to find epeak's error uses the same 
    propragation of error equation as any function.
    The matrix version of the equation is: 
    variance = PDMAT * COVARMAT * PDMAT.T
    uncertainty = np.sqrt(variance)
    Where PDMAT is the partial derivative matrix and COVARMAT is the 
    covariance matrix from the fit. 
    The function for computing epeak from ebreak or tem does not include 
    the amplitude term (i.e., norm__4 parameter), therefore there is no 
    apparent reason for including its partial derivative (which is 0). 
    However, we include it here so that the covariance matrix from the 
    model fit does not need to be changed. If we do edit the covariance 
    matrix and remove the norm__4 terms (which is row 3 and col 3 for 
    sbpl and grbm fits) we would get the same exact result as if we left 
    it in and had the partial deriv wrt norm as 0.

    """
    pars = None
    if model == 'grbm':
        # This looks confusing. The far left are the keys of the dict 
        #   and the right are the values. So 'alpha__1' is the key,
        #   'tem__3' is the value of the firt one. 
        #   The keys are the parameters and the values are their 
        #   respective function's partial deriv wrt those parameters. 
        pars = OrderedDict([
                    ('alpha__1', 
                        'tem__3'), 
                    ('beta__2', 
                        '0'), 
                    ('tem__3', 
                        'alpha__1 + 2.0'),
                    ('norm__4', 
                        '0')
                    ])
    elif model == 'cutoffpl':
        pars = OrderedDict([
                    ('PhoIndex__1', 
                        'HighECut__2'), 
                    ('HighECut__2', 
                        'PhoIndex__1 + 2.0'),
                    ('norm__3', 
                        '0')
                    ])
    elif model == 'sbpl':
        pars = OrderedDict([
            ('alpha__1',
                ('0.3*10**(0.3*atanh((alpha__1 + beta__2 + 4)/'
                 '(alpha__1 - beta__2)))*ebreak__3*(1/(alpha__1 '
                 '- beta__2) - (alpha__1 + beta__2 + 4)/(alpha__1 '
                 '- beta__2)**2)*log(10)/(1 - (alpha__1 + beta__2 '
                 '+ 4)**2/(alpha__1 - beta__2)**2)')
                ),
            ('beta__2',
                ('0.3*10**(0.3*atanh((alpha__1 + beta__2 + 4)/'
                 '(alpha__1 - beta__2)))*ebreak__3*(1/(alpha__1 - '
                 'beta__2) + (alpha__1 + beta__2 + 4)/(alpha__1 - '
                 'beta__2)**2)*log(10)/(1 - (alpha__1 + beta__2 + '
                 '4)**2/(alpha__1 - beta__2)**2)')
                ),
            ('ebreak__3',
                ('10**(0.3*atanh((alpha__1 + beta__2 + '
                 '4)/(alpha__1 - beta__2)))')
                ),
            ('norm__4', 
                        '0')
            ])
    else:
        raise Exception, "Model not recognized."
    return pars



