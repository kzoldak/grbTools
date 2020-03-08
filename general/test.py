
from grbTools.general.energetics import * #calc_flux, calc_fluence, calc_eiso


model = 'grbm+bbody+lpow'
parameter_dict_example(model)


from collections import OrderedDict
parameters = OrderedDict()

parameters['grbm'] = OrderedDict()
parameters['grbm']['alpha'] = -1.0
parameters['grbm']['beta'] = -2.5
parameters['grbm']['enterm'] = 500.0
parameters['grbm']['norm'] = 0.01
parameters['grbm']['entype'] = 'epeak' # optional
parameters['grbm']['enorm'] = 100.0    # optional

parameters['bbody'] = OrderedDict()
parameters['bbody']['kT'] = 52.8
parameters['bbody']['norm'] = 1.0
parameters['bbody']['program'] = 'xspec' # optional

parameters['lpow'] = OrderedDict()
parameters['lpow']['alpha'] = -1.2
parameters['lpow']['norm'] = 0.01
parameters['lpow']['enorm'] = 100.0    # optional

# parameters1 = dict({'beta': -2.5, 
#                             'alpha':-1.23, 
#                             'enterm': 306.9, 
#                             'norm': 0.0123})


flux = calc_flux(modelname=model, 
                 parameters=parameters,  
                 emin=10.0, 
                 emax=10000.0, 
                 redshift=None, 
                 photonflux=False)
print(flux)





# parameters1 = dict({'grbm':{'beta': -2.5, 
#                             'alpha':-1.23, 
#                             'enterm': 306.9, 
#                             'norm': 0.0123}})


# parameters2 = OrderedDict()
# parameters2['grbm'] = OrderedDict()
# parameters2['grbm']['alpha'] = -1.23
# parameters2['grbm']['beta'] = -2.5
# parameters2['grbm']['enterm'] = 306.9
# parameters2['grbm']['norm'] = 0.0123



# flux = calc_flux(modelname='grbm', 
#     parameters=parameters1,  
#     emin=10.0, emax=10000.0, redshift=None, 
#     photonflux=False)
# print(flux)


# flux = calc_flux(modelname='grbm', 
#     parameters=parameters2,  
#     emin=10.0, emax=10000.0, redshift=None, 
#     photonflux=False)
# print(flux)



