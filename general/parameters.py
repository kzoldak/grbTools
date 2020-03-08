
def parameter_dict_example(modelname):
    example = '''
    Notes
    -----
    - Anything labeled optional you do not HAVE to pass. The one shown is 
        the default. 
    - 'entype' options for 'grbm' and 'band' models: 'epeak', 'E0', 'tem', or 'ebreak'
    - 'entype' options for 'copl' model: 'epeak', 'E0', 'tem', 'epeak', or 'ebreak'
    - 'entype' options for 'sbpl' model: 'ebreak' or 'epeak'
    - 'grbm' and 'band' functions are identical, we proivde both for those who 
        prefer to use one name over the other. 
    - 'bbody' model has a 'program' attribute: 'xspec' or 'rmfit'.

    Parameter Dictionary Example for your model:
    --------------------------------------------

        from collections import OrderedDict
        parameters = OrderedDict()
    '''

    if 'grbm' in modelname:
        example = example + '''
        parameters['grbm'] = OrderedDict()
        parameters['grbm']['alpha'] = -1.0
        parameters['grbm']['beta'] = -2.5
        parameters['grbm']['enterm'] = 500.0
        parameters['grbm']['norm'] = 0.01
        parameters['grbm']['entype'] = 'epeak' # optional
        parameters['grbm']['enorm'] = 100.0    # optional
        '''
        
    if 'band' in modelname:
        example = example + '''
        parameters['band'] = OrderedDict()
        parameters['band']['alpha'] = -1.0
        parameters['band']['beta'] = -2.5
        parameters['band']['enterm'] = 500.0
        parameters['band']['norm'] = 0.01
        parameters['band']['entype'] = 'epeak' # optional
        parameters['band']['enorm'] = 100.0    # optional
        '''
    if 'sbpl' in modelname:
        example = example + '''
        parameters['sbpl'] = OrderedDict()
        parameters['sbpl']['alpha'] = -1.0
        parameters['sbpl']['beta'] = -2.5
        parameters['sbpl']['enterm'] = 500.0
        parameters['sbpl']['norm'] = 0.01
        parameters['sbpl']['entype'] = 'ebreak' # optional
        parameters['sbpl']['enorm'] = 100.0     # optional
        parameters['sbpl']['brkscale'] = 0.3    # optional
        '''
    if 'copl' in modelname:
        example = example + '''
        parameters['copl'] = OrderedDict()
        parameters['copl']['alpha'] = -1.0
        parameters['copl']['enterm'] = 500.0
        parameters['copl']['norm'] = 0.01
        parameters['copl']['entype'] = 'epeak' # optional
        parameters['copl']['enorm'] = 100.0    # optional
        '''
    if 'bbody' in modelname:
        example = example + '''
        parameters['bbody'] = OrderedDict()
        parameters['bbody']['kT'] = 52.8
        parameters['bbody']['norm'] = 1.0
        parameters['bbody']['program'] = 'xspec' # optional
        '''
    if 'lpow' in modelname:
        example = example + '''
        parameters['lpow'] = OrderedDict()
        parameters['lpow']['alpha'] = -1.2
        parameters['lpow']['norm'] = 0.01
        parameters['lpow']['enorm'] = 100.0    # optional
        '''
    if 'powerlaw' in modelname:
        example = example + '''
        parameters['powerlaw'] = OrderedDict()
        parameters['powerlaw']['alpha'] = 1.2
        parameters['powerlaw']['norm'] = 1.0
        parameters['powerlaw']['enorm'] = 1.0    # optional
        '''
    print(example)
