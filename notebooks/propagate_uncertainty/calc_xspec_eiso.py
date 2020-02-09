
from __future__ import division
from collections import OrderedDict
from grbTools.Energetics.eiso import Calc_Eiso

# cutoffpl+bbody v2  (GBM)
pars = OrderedDict()
pars['cutoffpl'] = OrderedDict()
pars['cutoffpl']['PhoIndex__1']= 1.2832349256355358
pars['cutoffpl']['HighECut__2']= 2151.6020054901182
pars['cutoffpl']['norm__3']= 4.058866228738239
pars['bbody'] = OrderedDict()
pars['bbody']['kT__4']= 47.36297544112143
pars['bbody']['norm__5']= 2.046087922228735


Eiso = Calc_Eiso(model='cutoffpl+bbody',
            parameters=pars, 
            emin = 10.,
            emax = 10000.,
            duration = 64.257-1.28,
            redshift = 4.35,
            program = 'xspec')
print(Eiso)




# grbm  (GBM)
pars = OrderedDict()
pars['grbm'] = OrderedDict()
pars['grbm']['alpha__1']= -1.0317748737257293
pars['grbm']['beta__2']= -2.2944235003425324
pars['grbm']['tem__3']= 536.8038210510679
pars['grbm']['norm__4']= 0.017453691740778158

Eiso = Calc_Eiso(model='grbm',
            parameters=pars, 
            emin = 10.,
            emax = 10000.,
            duration = 64.257-1.28,
            redshift = 4.35,
            program = 'xspec')
print(Eiso)



filename = '/Users/KimiZ/GRBs2/analysis/LAT/bn100728095/PYXSPEC/GBMwLAT/sbpl/xspec_fitresults_sbpl_-01-_L_.fit'

