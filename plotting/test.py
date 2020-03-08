import sys
sys.path.append('/Users/kimzoldak/GRBs/EisoCompare/')

import os

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
#%matplotlib inline

from matplotlib.backends.backend_pdf import PdfPages

from Zoldak.math.tools import root_mean_square as rms


from grbTools.plotting.journal import journal
journal()




f = ('/Users/kimzoldak/GRBs/sample/'
     'bestfit_model_vs_band__13_LAT_GRBs.txt')

df = pd.read_csv(f)


from grbTools.plotting.plot_comparisons import *

colnames = [col for col in df.columns if '_y' not in col]
dfa = df.loc[:, colnames].copy()

colnames = [col for col in df.columns if '_x' not in col]
dfb = df.loc[:, colnames].copy()

dfa.columns = dfa.columns.str.replace('_x', '')
dfb.columns = dfb.columns.str.replace('_y', '')



x = dfa['eiso'].apply(np.log10)
y = dfb['eiso'].apply(np.log10)

xerr = [dfa['eiso_ci_lo'].apply(np.log10), 
        dfa['eiso_ci_up'].apply(np.log10)]

yerr = [dfb['eiso_ci_lo'].apply(np.log10), 
        dfb['eiso_ci_up'].apply(np.log10)]

axlabel = r'$E_{iso}$ (erg)'
colors = ['blue',]
labels=['Band vs Best Model',]




# x-axis values (the Null Hypth) should always come second. 
deltas = np.asarray(dfb.eiso.apply(np.log10) - dfa.eiso.apply(np.log10))
delta = deltas.mean()
sigma = rms(deltas)


axLims = get_axes_limits(x, y, axBuffer=0.15)

#out_direc = '/Users/kimzoldak/Github/grbTools/plotting/'
out_direc = '/Users/kimzoldak/Documents/Thesis/Chapters/Figures/'
out_file = 'dummy.pdf'
out_file = os.path.join(out_direc, out_file)


pp = PdfPages(os.path.join(out_direc, out_file)) # must end in .pdf
plot_diagonal(axLims)
plot_comparison(x=x, y=x, xerr=xerr, yerr=yerr, parameter='eiso', xaxlabel='Band', 
    yaxlabel='Best Model', axLims=axLims,)
plot_stats(n=len(x), delta=delta, sigma=sigma, xloc=0.17, yloc=0.85)   # 0.05, 0.93
#plt.subplots_adjust(left=0.15, right=0.96, top=0.96, bottom=0.14)
plt.savefig(pp, format='pdf', bbox_inches='tight')
plt.close()
pp.close()

os.system('open %s'%out_file)




