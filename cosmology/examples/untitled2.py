#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 12:33:13 2020

@author: kimzoldak
"""

import numpy as np
from numpy import pi

import pandas as pd

import astropy
from astropy import units as u

import matplotlib.pyplot as plt
%matplotlib inline


#def get_moduli(model, redshifts=None, H0=None, Om=None):
#    """
#    model : str. 'concord', 'weyl', or 'riess'
#                    or 'c', 'w', and 'r'
#   
#    """
#    H0 = H_knot
#    OM = omega_m
#    if model.startswith('c'):  
#        DLs = [LumDist_concordance(redshift=z, H_knot=H0, omega_m=OM) for z in redshifts]
#    elif model.startswith('r'):  
#        DLs = [LumDist_riess(redshift=z, H_knot=H0, omega_m=OM) for z in redshifts]
#    elif model.startswith('w'):  
#        DLs = [LumDist_weylgravity(redshift=z, H_knot=H0) for z in redshifts]
#    moduli = [distance_modulus(lumdistance=dL) for dL in DLs]
#    return moduli


def distance_modulus(lumdistances):
    dls = np.asarray(lumdistances)
    return 5*np.log10(dls)-5

def calc_eiso(fluences, redshifts, lumdistances):
    S = np.asarray(fluences)
    DL = np.asarray(lumdistances)
    z = np.asarray(redshifts)
    return S*((4.0*pi*(DL**2))/(1.0+z))



#```
#Subscripts:
#1) LCDM:          H0 = 67.8, Om = 0.308
#2) Riess:         H0 = 67.8, Om = 0.308
#3) Weyl Gravity:  H0 = 67.8
#4) LCDM:          H0 = 65.0, Om = 0.308
#5) LCDM:          H0 = 75.0, Om = 0.308
#```
    
df1 = pd.read_csv('eisoenergies_diff_cosmo.txt', sep=',', skiprows=1)
df2 = pd.read_csv('eisoenergies_diff_H0.txt', sep=',', skiprows=1)
# remove duplicate to eiso1 in our other dataframe.
[df2.pop(col) for col in df2.columns if 'eiso1' in col];
df2.pop('DL1');
# rename columns 
colnames = df2.columns
colnames = colnames.str.replace('eiso2', 'eiso4')
colnames = colnames.str.replace('eiso3', 'eiso5')
colnames = colnames.str.replace('DL2', 'DL4')
colnames = colnames.str.replace('DL3', 'DL5')
df2.columns = colnames
# merge dataframes and keep only unique columns. The rest are repeats,
# hence our need to rename the cols we wanted to keep. 
df = pd.merge(df1, df2)



H0 = 67.8
Om = 0.308

distance_modulus(lumdistances)

# Concordance Cosmology (Lambda CDM)
df['eiso1']         = calc_eiso(fluence=df.fluence, lumdist=df['DL1'], redshift=df.z)
df['eiso1_err_low'] = calc_eiso(fluence=df.fluence_err_low, lumdist=df['DL1'], redshift=df.z)
df['eiso1_err_up']  = calc_eiso(fluence=df.fluence_err_up, lumdist=df['DL1'], redshift=df.z)
df['eiso1_err']     = df['eiso1']-df['eiso1_err_low']  # margin of error. Eiso errs should be symmetrical when unlogged

# Riess Cosmology
df['eiso2']         = calc_eiso(fluence=df.fluence, lumdist=df['DL2'], redshift=df.z)
df['eiso2_err_low'] = calc_eiso(fluence=df.fluence_err_low, lumdist=df['DL2'], redshift=df.z)
df['eiso2_err_up']  = calc_eiso(fluence=df.fluence_err_up, lumdist=df['DL2'], redshift=df.z)
df['eiso2_err']     = df['eiso2']-df['eiso2_err_low']  # margin of error. Eiso errs should be symmetrical when unlogged


# Weyl Gravity
df['eiso3']         = calc_eiso(fluence=df.fluence, lumdist=df['DL3'], redshift=df.z)
df['eiso3_err_low'] = calc_eiso(fluence=df.fluence_err_low, lumdist=df['DL3'], redshift=df.z)
df['eiso3_err_up']  = calc_eiso(fluence=df.fluence_err_up, lumdist=df['DL3'], redshift=df.z)
df['eiso3_err']     = df['eiso3']-df['eiso3_err_low']  # margin of error. Eiso errs should be symmetrical when unlogged



axLims = (51.51, 55.2)


plt.clf()
plt.figure(figsize=(8,7))
#plt.figure()


xaxis = 'eiso1'
yaxis = 'eiso2'
pltKwgs = dict(fmt='o', color='white', ecolor='blue', ms=5, lw=0.6, mec='blue', mew=0.75,
               capsize=0, alpha=1, label='x: LCDM\ny: Riess')
x = df[xaxis].apply(np.log10)
y = df[yaxis].apply(np.log10)
xerr = np.asarray([(df[xaxis].apply(np.log10)-df[xaxis+'_err_low'].apply(np.log10)).values, # lower margin of error
                   (df[xaxis+'_err_up'].apply(np.log10)-df[xaxis].apply(np.log10)).values]) # upper margin of error
yerr = np.asarray([(df[yaxis].apply(np.log10)-df[yaxis+'_err_low'].apply(np.log10)).values, # lower margin of error
                   (df[yaxis+'_err_up'].apply(np.log10)-df[yaxis].apply(np.log10)).values]) # upper margin of error
#plt.errorbar(x=x, y=y, yerr=yerr, xerr=xerr, **pltKwgs)
markers, caps, bars = plt.errorbar(x=x, y=y, yerr=yerr, xerr=xerr, **pltKwgs)
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]


xaxis = 'eiso1'
yaxis = 'eiso3'
pltKwgs = dict(fmt='o', color='white', ecolor='red', ms=5, lw=0.6, mec='red', mew=0.75,
               capsize=0, alpha=1, label='x: LCDM\ny: Weyl')
x = df[xaxis].apply(np.log10)
y = df[yaxis].apply(np.log10)
xerr = np.asarray([(df[xaxis].apply(np.log10)-df[xaxis+'_err_low'].apply(np.log10)).values, # lower margin of error
                   (df[xaxis+'_err_up'].apply(np.log10)-df[xaxis].apply(np.log10)).values]) # upper margin of error
yerr = np.asarray([(df[yaxis].apply(np.log10)-df[yaxis+'_err_low'].apply(np.log10)).values, # lower margin of error
                   (df[yaxis+'_err_up'].apply(np.log10)-df[yaxis].apply(np.log10)).values]) # upper margin of error
#plt.errorbar(x=x, y=y, yerr=yerr, xerr=xerr, **pltKwgs)
markers, caps, bars = plt.errorbar(x=x, y=y, yerr=yerr, xerr=xerr, **pltKwgs)
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]
   

xaxis = 'eiso2'
yaxis = 'eiso3'
pltKwgs = dict(fmt='o', color='white', ecolor='green', ms=5, lw=0.6, mec='green', mew=0.75,
               capsize=0, alpha=1, label='x: Riess\ny: Weyl')
x = df[xaxis].apply(np.log10)
y = df[yaxis].apply(np.log10)
xerr = np.asarray([(df[xaxis].apply(np.log10)-df[xaxis+'_err_low'].apply(np.log10)).values, # lower margin of error
                   (df[xaxis+'_err_up'].apply(np.log10)-df[xaxis].apply(np.log10)).values]) # upper margin of error
yerr = np.asarray([(df[yaxis].apply(np.log10)-df[yaxis+'_err_low'].apply(np.log10)).values, # lower margin of error
                   (df[yaxis+'_err_up'].apply(np.log10)-df[yaxis].apply(np.log10)).values]) # upper margin of error
#plt.errorbar(x=x, y=y, yerr=yerr, xerr=xerr, **pltKwgs)
markers, caps, bars = plt.errorbar(x=x, y=y, yerr=yerr, xerr=xerr, **pltKwgs)
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]
#[marker.set_alpha(0.25) for marker in markers]


# Plot a line at equal values of x and y
xlinedata = np.linspace(axLims[0], axLims[1], 10)
ylinedata = 1.0*xlinedata+0
plt.plot(xlinedata, ylinedata, 'k-', lw=1, alpha=0.25)

plt.xlim(*axLims)
plt.ylim(*axLims)
plt.legend(loc=0, fontsize=11, labelspacing=0.7, handletextpad=0, frameon=False)
plt.xlabel('$E_{iso}$ (erg)', fontsize=13)
plt.ylabel('$E_{iso}$ (erg)', fontsize=13)
plt.show()



df['diff_1v2'] = (df['eiso1'].apply(np.log10) - df['eiso2'].apply(np.log10)).apply(abs)

df['diff_1v3'] = (df['eiso1'].apply(np.log10) - df['eiso3'].apply(np.log10)).apply(abs)

df['diff_2v3'] = (df['eiso2'].apply(np.log10) - df['eiso3'].apply(np.log10)).apply(abs)

sort_by = 'diff_1v2'
kim = df.sort_values(by=sort_by).loc[:, ['z',sort_by]]

plt.clf()
(kim).plot(kind='scatter', x='z', y=sort_by, figsize=(8,5),
                                                        marker='o', s=50, color='white', edgecolor='blue')


df.sort_values(by='z').loc[:, ['z', 'number', 'LATburst']]


sort_by = 'diff_1v3'
kim2 = df.sort_values(by=sort_by).loc[:, ['z',sort_by]]

sort_by = 'diff_2v3'
kim3 = df.sort_values(by=sort_by).loc[:, ['z',sort_by]]




xLims = (None, None)
yLims = (None, None)

# xLims = (0.2, 1)
# yLims = (-0.001, .05)


plt.clf()
plt.figure(figsize=(8,7))
#plt.figure()

plt.plot(df['z'], df['diff_1v2'], marker='o', ms=8, color='white', mec='blue', mew=2, lw=0,  label='x: LCDM\ny: Riess')

plt.plot(df['z'], df['diff_1v3'], marker='o', ms=8, color='white', mec='red', mew=2, lw=0, label='x: LCDM\ny: Weyl')

plt.plot(df['z'], df['diff_2v3'], marker='o', ms=8, color='white', mec='green', mew=2, lw=0, label='x: Riess\ny: Weyl')


plt.xlim(*xLims)
plt.ylim(*yLims)
plt.legend(loc=0, fontsize=11, labelspacing=0.7, handletextpad=0, frameon=False)
plt.ylabel('Diff $\log_{10}$($E_{iso}$) (erg)', fontsize=13)
plt.xlabel('$z$', fontsize=13)
#plt.yscale('log')
#plt.xscale('log')
plt.show()
