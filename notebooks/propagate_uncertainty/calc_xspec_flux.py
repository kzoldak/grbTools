'''
This program shows you how to calculate energy flux from an xspec model 
fit. We build a dictionary of the model and parameters and then pass it 
to Calc_Flux. This particular Calc_Flux function was set up for XSPEC 
models only. 

pars dictionary
---------------
pars is an ordered dictionary, meaning the keys and values stay in the 
  order they were entered. This is critical for this program because the 
  Calc_Flux function was set up so that a list of parameters is passed 
  to the model functions, not a dictionary. 

  For example:
  grbm(energy, alpha, beta, tem, norm) 
  is the model set up for the band function. 

  lpow(plIndex, norm) 
  is the model set up for the power-law function. 

  Within the pars dictionary, the parameters must be set up as so:
      pars = OrderedDict()

      pars['grbm'] = OrderedDict()
      pars['grbm']['alpha__1']= -1.0
      pars['grbm']['beta__2']= -2.3
      pars['grbm']['tem__3']= 300.0
      pars['grbm']['norm__4']= 0.01

      pars['lpow'] = OrderedDict()
      pars['lpow']['plIndex__5']= -1.0
      pars['lpow']['norm__6']= 0.001

  Notice the order of the parameters. 
  If we do 
      pars['grbm']['beta__2']= -2.3
      pars['grbm']['alpha__1']= -1.0
      pars['grbm']['tem__3']= 300.0
      pars['grbm']['norm__4']= 0.01
  instead, then the band function reads the -2.3 as alpha, not beta. 
  Thus, having the parameter names defined within the dictionary is 
  not there to help the model functions (being integrated) know the 
  proper parameter values, but to keep things organized on our end. 
  
  You could do the following, and it wouldn't matter:
      pars['grbm']['alp']= -2.3
      pars['grbm']['bet']= -1.0
      pars['grbm']['ech']= 300.0
      pars['grbm']['amp']= 0.01
  But we don't recommend doing this. 


  Setting up models
  ------------------
  Regardless of whether or not you have additive models, you should 
  always define the model name as the first level of keys in the nested 
  dictionary. 
  For example, if we only have the band function results, we would still 
  do:
        pars['grbm'] = OrderedDict()
        pars['grbm']['alpha__1'] = -1.0
        pars['grbm']['beta__2'] = -2.3
        pars['grbm']['tem__3'] = 300.0
        pars['grbm']['norm__4'] = 0.01

'''
from __future__ import division
from collections import OrderedDict
from grbTools.XSPEC.flux import Calc_Flux

#from grbTools.Energetics.flux import Calc_Flux


burst = 'bn080916009'
z = 4.35


# grbm+lpow
pars = OrderedDict()
pars['grbm'] = OrderedDict()
pars['grbm']['alpha__1']= -0.7429962860863923
pars['grbm']['beta__2']= -2.177873923426739
pars['grbm']['tem__3']= 305.50532574626766
pars['grbm']['norm__4']= 0.019014347257507447
pars['lpow'] = OrderedDict()
pars['lpow']['plIndex__5']= -2.0771292778765464
pars['lpow']['norm__6']= 0.0009925643425563427

Flux = Calc_Flux(model='grbm+lpow', 
            parameters=pars, 
            emin=10.0, 
            emax=10000.0, 
            redshift=4.35, 
            photonflux=False)
# 1.4642929905655007e-06
print(Flux)


# sbpl+bbody v2
pars = OrderedDict()
pars['sbpl'] = OrderedDict()
pars['sbpl']['alpha__1']= -1.2859259039727824
pars['sbpl']['beta__2']= -2.230432546708862
pars['sbpl']['ebreak__3']= 574.7141875888526
pars['sbpl']['norm__4']= 0.010828176802794976
pars['bbody'] = OrderedDict()
pars['bbody']['kT__5']= 45.39075663906129
pars['bbody']['norm__6']= 1.7800897291179052

Flux = Calc_Flux(model='sbpl+bbody', 
            parameters=pars, 
            emin=10.0, 
            emax=10000.0, 
            redshift=4.35, 
            photonflux=False)
# 1.521006779408721e-06
print(Flux)


# cutoffpl+bbody v2  (GBM)
pars = OrderedDict()
pars['cutoffpl'] = OrderedDict()
pars['cutoffpl']['PhoIndex__1']= 1.2832349256355358
pars['cutoffpl']['HighECut__2']= 2151.6020054901182
pars['cutoffpl']['norm__3']= 4.058866228738239
pars['bbody'] = OrderedDict()
pars['bbody']['kT__4']= 47.36297544112143
pars['bbody']['norm__5']= 2.046087922228735

Flux = Calc_Flux(model='cutoffpl+bbody', 
            parameters=pars, 
            emin=10.0, 
            emax=10000.0, 
            redshift=4.35, 
            photonflux=False)
# 1.600751338927506e-06
print(Flux)



# pars = OrderedDict()
# pars['cutoffpl'] = OrderedDict()
# pars['cutoffpl']['PhoIndex__1']= 1.2832349256355358
# pars['cutoffpl']['HighECut__2']= 2151.6020054901182
# pars['cutoffpl']['norm__3']= 4.058866228738239
# pars['bbody'] = OrderedDict()
# pars['bbody']['kT__4']= 47.36297544112143
# pars['bbody']['norm__5']= 2.046087922228735

# Flux = Calc_Flux(model='cutoffpl+bbody', 
#             parameterdict=pars, 
#             emin=10.0, 
#             emax=10000.0, 
#             redshift=4.35, 
#             program='xspec',
#             photonflux=False)
# print(Flux)




