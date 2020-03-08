import sys
#sys.path.append('/Users/kimzoldak/GRBs/EisoCompare/')

import os

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

from grbTools.math.tools import root_mean_square as rms

from grbTools.compare.plot_comparisons import *  # imports journal for us.


df = pd.read_csv(('/Users/kimzoldak/GRBs/sample/'
                    'band_sbpl_28_overlap.txt'), skiprows=4)


# x-axis  (dfa)
colnames = [col for col in df.columns if '_y' not in col]
dfa = df.loc[:, colnames].copy()

# y-axis  (dfb)
colnames = [col for col in df.columns if '_x' not in col]
dfb = df.loc[:, colnames].copy()


dfa.columns = dfa.columns.str.replace('_x', '')
dfb.columns = dfb.columns.str.replace('_y', '')



pltKwargs = [dict(fmt = 'o', ms=0, capsize = None, color = 'white', 
                  ecolor = 'white', mec='white', elinewidth=0, alpha=0.0),
             
             dict(fmt = 'o', ms=5, capsize = None, color = 'blueviolet', 
                  ecolor = 'blueviolet', mec='k', elinewidth=1, alpha=0.7),
             
             dict(fmt = 'o', ms=5, capsize = None, color = 'plum', 
                  ecolor = 'plum', mec='mediumorchid', elinewidth=1, alpha=0.5)]




out_file = 'dummy.pdf'
pp = PdfPages(out_file)


for par in ['eiso','epeak','epeakRest','alpha', 'beta', 'norm', 'flux']:
    plt.clf()
    if par in ['alpha', 'beta']:
        func = lambda XVAR: XVAR*1.0  # leave linear
    else:
        func = np.log10       # log the data

    x = dfa[par].apply(func)
    y = dfb[par].apply(func)
    xerr = [dfa[par+'_ci_lo'].apply(func), 
            dfa[par+'_ci_up'].apply(func)]
    yerr = [dfb[par+'_ci_lo'].apply(func), 
            dfb[par+'_ci_up'].apply(func)]

    # x-axis values (the Null Hypth) should always come second. 
    deltas = np.asarray(dfb[par].apply(func) - dfa[par].apply(func))
    delta = deltas.mean()
    sigma = rms(deltas)

    axLims = get_axes_limits(x, y, axBuffer=0.15)

    plot_diagonal(axLims)
    
    plot_comparison(x=x, y=y, xerr=xerr, yerr=yerr, parameter=par, 
        xaxlabel='Band', 
        yaxlabel='Sbpl', 
        axLims=axLims, axBuffer=0.15, ax=None, 
        **pltKwargs[1])
    
    plot_stats(n=len(x), delta=delta, sigma=sigma, xloc=0.17, yloc=0.85)   # 0.05, 0.93
    
    plt.subplots_adjust(left=0.15, right=0.96, top=0.96, bottom=0.14)
    plt.show()
    #plt.savefig(pp, format='pdf', )#bbox_inches='tight')
    print(par)

    #plt.close()
# pp.close()

# os.system('open %s'%out_file)


# os.system('open ./')





# for par in ['eiso','epeak','epeakRest','alpha', 'beta', 'norm', 'flux']:
#     if par in ['alpha', 'beta']:
#         func = 1  # leave linear
#     else:
#         func = np.log10  # log the data

#         x = dfa[par]
#         y = dfb[par]
#         xerr = [dfa[par+'_ci_lo'], 
#                 dfa[par+'_ci_up']]
#         yerr = [dfb[par+'_ci_lo'], 
#                 dfb[par+'_ci_up']]
#     else:
#         x = dfa[par].apply(np.log10)
#         y = dfb[par].apply(np.log10)
#         xerr = [dfa['eiso_ci_lo'].apply(np.log10), 
#                 dfa['eiso_ci_up'].apply(np.log10)]
#         yerr = [dfb['eiso_ci_lo'].apply(np.log10), 
#                 dfb['eiso_ci_up'].apply(np.log10)]


