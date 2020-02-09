#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 13:07:44 2020

@author: kimzoldak
"""

import numpy as np
from numpy import pi

import pandas as pd

import astropy
from astropy import units as u

import matplotlib.pyplot as plt
%matplotlib inline

from lumdist import Cosmology as Cosmo
cosmo = Cosmo()


def distance_modulus(lumdistances):
    dls = np.asarray(lumdistances)
    return 5*np.log10(dls)-5


def distance_modulus_2(lumdistances):
    dls = np.asarray(lumdistances)
    return 25 + 5*np.log10(dls)



def calc_eiso(fluences, redshifts, lumdistances):
    S = np.asarray(fluences)
    DL = np.asarray(lumdistances)
    z = np.asarray(redshifts)
    return S*((4.0*pi*(DL**2))/(1.0+z))

deltas = lambda x,y: np.asarray(abs(np.log10(x)-np.log10(y)))


def plot_hubble_diagram(x, y, xLims=None, yLims=None, ax=None, **pltKwgs):
    # plt.rcParams 
    # plt.rcParams['ytick.right'] = False

    if plt.rcParams['ytick.right'] is False:
        plt.rcParams['ytick.right'] = True
    
    if ax is None:
        ax = plt.gca()
    if not pltKwgs:
        pltKwgs = dict(ls='--', color='blue', alpha=0.5)
    PLT = ax.plot(x, y, **pltKwgs)
    ax.minorticks_on()
    ax.set_xlim(*xLims)
    ax.set_ylim(*yLims)
    ax.set_xlabel(r'$z$', fontsize=13)
    ax.set_ylabel(r'$\mu$ (mag)', fontsize=13)
    ax.legend(loc=0, ) #fontsize=11, labelspacing=0.7, handletextpad=0, frameon=False)
    plt.rcParams['ytick.right'] = False
    return PLT




# General Hubble Diagram with our models.
redshifts = np.linspace(0.1, 8, 100)

# FlatLambdaCDM 
cosmoconstants = {'H0':67.8, 'Om': 0.308}
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='FlatLambdaCDM', 
                                 units='pc',
                                 **cosmoconstants)
mag1 = distance_modulus(lumdistances=DL)
del DL


# w0wzCDM   with 'w0': -1.31, 'wp': 1.48
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-1.31, 'wp':1.48}
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0wpCDM', 
                                 units='pc',
                                 **cosmoconstants)
mag2 = distance_modulus(lumdistances=DL)

# w0waCDM   with 'w0':-0.9, 'wa':0.2

# w0=0.2, wa=-1.4
# w0=0.1, wa=-0.9
#cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-0.9, 'wa':0.2} 
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':0.2, 'wa':-0.9}    
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0waCDM', 
                                 units='pc',
                                 **cosmoconstants)
mag3 = distance_modulus(lumdistances=DL)
del DL


# weylgravity
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'q0':-0.37}
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='weylgravity', 
                                 units='pc',
                                 **cosmoconstants)
mag4 = distance_modulus(lumdistances=DL)
del DL


# Hubble Diagram
colors = ['blue', 'red', 'green', 'orange']
#labels = ['mag1', 'mag2', 'mag3', 'mag4']
labels = ['FlatLCDM', 'w0wp', 'w0wa', 'WeylGravity']


plt.clf()
plt.figure(figsize=(6, 4.5))
for mag,clr,lab in zip([mag1, mag2, mag3, mag4], colors, labels):
    pltKwgs = dict(ls='-', alpha=0.75, lw=2)
    pltKwgs['color'] = clr
    pltKwgs['label'] = lab
    plot_hubble_diagram(redshifts, mag, 
                        #xLims=(0, 5.5),
                        xLims=(0, 7), 
                        yLims=(35, 50), 
                        ax=None, 
                        **pltKwgs)
#plt.show()
#plt.savefig('3.pdf')
    
    
wa = -0.9
w0 = 0.1
wa + 2.2*w0 , wa-0.45*w0



plt.clf()
plt.figure(figsize=(10, 9))
for mag,clr,lab in zip([mag1, mag2, mag3, mag4], colors, labels):
    pltKwgs = dict(ls='-', alpha=0.75, lw=2)
    pltKwgs['color'] = clr
    pltKwgs['label'] = lab
    plot_hubble_diagram(redshifts, mag, 
                        #xLims=(0, 5.5),
                        xLims=(0, 1), 
                        yLims=(38, 44), 
                        ax=None, 
                        **pltKwgs)
    
plt.show()

#plt.savefig('3.pdf')

#
#def plot_hubble_diagram(yaxes, colors, labels, xLims=None, yLims=None, 
#                        ax=None, **pltKwgs):
#    if ax is None:
#        ax = plt.gca()
#    
#    # Plot a line at equal values of x and y
#    #xlinedata = np.linspace(axLims[0], axLims[1], 10)
#    #ylinedata = 1.0*xlinedata+0
#    #PLT = ax.plot(xlinedata, ylinedata, 'k-', lw=1, alpha=0.25)
#    x = df.z
#    for yaxis,color,label in zip(yaxes,colors,labels):
#        if not pltKwgs:
#            pltKwgs = dict(ls='--', color='blue', alpha=0.5)
#        pltKwgs['mec'] = color
#        pltKwgs['label'] = label
#        y = df[yaxis]
#        PLT = ax.plot(x, y, **pltKwgs)
#    ax.set_xlim(*xLims)
#    ax.set_ylim(*yLims)
#    ax.set_xlabel('$z$', fontsize=13)
#    ax.set_ylabel(r'$\mu$ (mag)', fontsize=13)
#    ax.legend(loc=0, ) #fontsize=11, labelspacing=0.7, handletextpad=0, frameon=False)
#    return PLT

    
def plot_differences(xaxes, yaxes, axLims, colors, labels, axlabel, ax=None):
    if ax is None:
        ax = plt.gca()
    
    # Plot a line at equal values of x and y
    xlinedata = np.linspace(axLims[0], axLims[1], 10)
    ylinedata = 1.0*xlinedata+0
    PLT = ax.plot(xlinedata, ylinedata, 'k-', lw=1, alpha=0.25)
    for xaxis,yaxis,color,label in zip(xaxes,yaxes,colors,labels):
        pltKwgs = dict(marker='o', color='white', ms=5, lw=0.6, mew=0.75, alpha=1)
        pltKwgs['mec'] = color
        pltKwgs['label'] = label
        x = df[xaxis].apply(np.log10)
        y = df[yaxis].apply(np.log10)
        ax.plot(x, y, **pltKwgs)
    ax.set_xlim(*axLims)
    ax.set_ylim(*axLims)
    ax.set_xlabel(axlabel, fontsize=13)
    ax.set_ylabel(axlabel, fontsize=13)
    ax.legend(loc=0, ) #fontsize=11, labelspacing=0.7, handletextpad=0, frameon=False)
    return PLT





# Hubble Diagram
#xaxes = ['eiso1', 'eiso1', 'eiso1']
yaxes = ['mag1', 'mag2', 'mag3', 'mag4']
colors = ['blue', 'red', 'green', 'orange']
labels = ['a','b','c', 'd']

plt.clf()
plot_hubble_diagram(yaxes=yaxes, 
                    colors=colors, 
                    labels=labels, 
                    xLims = (0, 5.5),
                    yLims = (35, 50),
                    ax=None)
plt.show()







df = pd.read_csv('eisoenergies_diff_cosmo.txt', sep=',', skiprows=1)

df = df.loc[:, :'fluence_err_up']


#```
#Subscripts:
#1) LCDM:          H0 = 67.8, Om = 0.308
#2) Riess:         H0 = 67.8, Om = 0.308
#3) Weyl Gravity:  H0 = 67.8
#4) LCDM:          H0 = 65.0, Om = 0.308
#5) LCDM:          H0 = 75.0, Om = 0.308
#```
    

redshifts = np.asarray(df.z)


# FlatLambdaCDM 
cosmoconstants = {'H0':67.8, 'Om': 0.308}
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='FlatLambdaCDM', 
                                 units='cm',
                                 **cosmoconstants)
df['DL1'] = DL
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='FlatLambdaCDM', 
                                 units='pc',
                                 **cosmoconstants)
df['mag1'] = distance_modulus(lumdistances=DL)
zone1 = cosmo.luminosity_distance(redshifts=[1.0,], 
                                 model='FlatLambdaCDM', 
                                 units='cm',
                                 **cosmoconstants)
del DL


# w0wzCDM   with 'w0': -1.31, 'wp': 1.48
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-1.31, 'wp':1.48}
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0wpCDM', 
                                 units='cm',
                                 **cosmoconstants)
df['DL2'] = DL
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0wpCDM', 
                                 units='pc',
                                 **cosmoconstants)
df['mag2'] = distance_modulus(lumdistances=DL)
zone2 = cosmo.luminosity_distance(redshifts=[1.0,], 
                                 model='w0wpCDM', 
                                 units='cm',
                                 **cosmoconstants)



# w0waCDM   with 'w0':-0.9, 'wa':0.2
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-0.9, 'wa':0.2}
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0waCDM', 
                                 units='cm',
                                 **cosmoconstants)
df['DL3'] = DL
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0waCDM', 
                                 units='pc',
                                 **cosmoconstants)
df['mag3'] = distance_modulus(lumdistances=DL)
zone3 = cosmo.luminosity_distance(redshifts=[1.0,], 
                                 model='w0waCDM', 
                                 units='cm',
                                 **cosmoconstants)
del DL


# weylgravity
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'q0':-0.37}
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='weylgravity', 
                                 units='cm',
                                 **cosmoconstants)
df['DL4'] = DL
DL = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='weylgravity', 
                                 units='pc',
                                 **cosmoconstants)
df['mag4'] = distance_modulus(lumdistances=DL)
zone4 = cosmo.luminosity_distance(redshifts=[1.0,], 
                                 model='weylgravity', 
                                 units='cm',
                                 **cosmoconstants)
del DL





#EISO1 = calc_eiso(fluences=df.fluence, redshifts=df.z, lumdistances=df.DL1)
#EISO2 = calc_eiso(fluences=df.fluence, redshifts=df.z, lumdistances=df.DL2)
#EISO3 = calc_eiso(fluences=df.fluence, redshifts=df.z, lumdistances=df.DL3)
#EISO4 = calc_eiso(fluences=df.fluence, redshifts=df.z, lumdistances=df.DL4)
df['eiso1'] = calc_eiso(fluences=df.fluence, redshifts=df.z, lumdistances=df.DL1)
df['eiso2'] = calc_eiso(fluences=df.fluence, redshifts=df.z, lumdistances=df.DL2)
df['eiso3'] = calc_eiso(fluences=df.fluence, redshifts=df.z, lumdistances=df.DL3)
df['eiso4'] = calc_eiso(fluences=df.fluence, redshifts=df.z, lumdistances=df.DL4)



# Differences
df['Diff_DL_1v2'] = deltas(df.DL1, df.DL2)
df['Diff_DL_1v3'] = deltas(df.DL1, df.DL3)
df['Diff_DL_1v4'] = deltas(df.DL1, df.DL4)


df['Diff_eiso_1v2'] = deltas(df.eiso1, df.eiso2)
df['Diff_eiso_1v3'] = deltas(df.eiso1, df.eiso3)
df['Diff_eiso_1v4'] = deltas(df.eiso1, df.eiso4)

df['Diff_eiso_2v1'] = deltas(df.eiso2, df.eiso1)
df['Diff_eiso_2v3'] = deltas(df.eiso2, df.eiso3)
df['Diff_eiso_2v4'] = deltas(df.eiso2, df.eiso4)

df['Diff_eiso_3v1'] = deltas(df.eiso3, df.eiso1)
df['Diff_eiso_3v2'] = deltas(df.eiso3, df.eiso2)
df['Diff_eiso_3v4'] = deltas(df.eiso3, df.eiso4)

df['Diff_eiso_4v1'] = deltas(df.eiso4, df.eiso1)
df['Diff_eiso_4v2'] = deltas(df.eiso4, df.eiso2)
df['Diff_eiso_4v3'] = deltas(df.eiso4, df.eiso3)




# Hubble Diagram
#xaxes = ['eiso1', 'eiso1', 'eiso1']
yaxes = ['mag1', 'mag2', 'mag3', 'mag4']
colors = ['blue', 'red', 'green', 'orange']
labels = ['a','b','c', 'd']

plt.clf()
plot_hubble_diagram(yaxes=yaxes, 
                    colors=colors, 
                    labels=labels, 
                    xLims = (0, 5.5),
                    yLims = (35, 50),
                    ax=None)
plt.show()




xaxes = ['eiso1', 'eiso1', 'eiso1']
yaxes = ['eiso2', 'eiso3', 'eiso4']
colors = ['blue', 'red', 'green']
labels = ['x: LCDM\ny: w0wp', 'x: LCDM\ny: w0wa', 'x: LCDM\ny: Weyl']
axlabel = r'$E_{iso}$ (erg)'

plt.clf()
plot_differences(xaxes=xaxes, yaxes=yaxes, axLims=(51.51, 55.2), 
                 colors=colors, 
                 labels=labels, 
                 axlabel=axlabel,
                 ax=None)
plt.show()

# 1) LCDM
# 2) w0wp
# 3) w0wa
# 4) weylgravity

xaxes = ['eiso1', 'eiso1', 'eiso1', ]
yaxes = ['eiso2', 'eiso3', 'eiso4']
colors = ['blue', 'red', 'green']
labels = ['x: LCDM\ny: w0wp', 'x: LCDM\ny: w0wa', 'x: LCDM\ny: Weyl']
axlabel = r'$E_{iso}$ (erg)'

plt.clf()
plot_differences(xaxes=xaxes, yaxes=yaxes, axLims=(51.51, 55.2), 
                 colors=colors, 
                 labels=labels, 
                 axlabel=axlabel,
                 ax=None)
plt.show()








sort_by = 'Diff_Eiso_1v2'
kim = df.sort_values(by=sort_by).loc[:, ['z',sort_by]]
plt.clf()
(kim).plot(kind='scatter', x='z', y=sort_by, 
        figsize=(8,5), marker='o', s=50, color='white', edgecolor='blue')


sort_by = 'Diff_Eiso_1v3'
kim2 = df.sort_values(by=sort_by).loc[:, ['z',sort_by]]
plt.clf()
(kim2).plot(kind='scatter', x='z', y=sort_by, 
           figsize=(8,5), marker='o', s=50, color='white', edgecolor='blue')



sort_by = 'Diff_Eiso_1v4'
kim3 = df.sort_values(by=sort_by).loc[:, ['z',sort_by]]
plt.clf()
(kim3).plot(kind='scatter', x='z', y=sort_by, 
           figsize=(8,5), marker='o', s=50, color='white', edgecolor='blue')






axLims = (51.51, 55.2)

plt.clf()
plt.figure(figsize=(8,7))
#plt.figure()

xaxis = 'eiso1'
yaxis = 'eiso2'
pltKwgs = dict(marker='o', color='white', ms=5, lw=0.6, mew=0.75, mec='blue',
               alpha=1, label='x: LCDM\ny: w0wp')
x = df[xaxis].apply(np.log10)
y = df[yaxis].apply(np.log10)
plt.plot(x, y, **pltKwgs)

xaxis = 'eiso1'
yaxis = 'eiso3'
pltKwgs = dict(marker='o', color='white', ms=5, lw=0.6, mew=0.75, mec='red',
               alpha=1, label='x: LCDM\ny: w0wa')
x = df[xaxis].apply(np.log10)
y = df[yaxis].apply(np.log10)
plt.plot(x, y, **pltKwgs)
   

xaxis = 'eiso1'
yaxis = 'eiso4'
pltKwgs = dict(marker='o', color='white', ms=5, lw=0.6, mew=0.75, mec='green',
               alpha=1, label='x: LCDM\ny: Weyl')
x = df[xaxis].apply(np.log10)
y = df[yaxis].apply(np.log10)
plt.plot(x, y, **pltKwgs)

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









axLims = (51.51, 55.2)


plt.clf()
plt.figure(figsize=(8,7))
#plt.figure()

xaxis = 'eiso1'
yaxis = 'eiso2'
pltKwgs = dict(fmt='o', color='white', ecolor='blue', ms=5, lw=0.6, mec='blue', mew=0.75,
               capsize=0, alpha=1, label='x: LCDM\ny: w0wp')
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
               capsize=0, alpha=1, label='x: LCDM\ny: w0wa')
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
yaxis = 'eiso4'
pltKwgs = dict(fmt='o', color='white', ecolor='green', ms=5, lw=0.6, mec='green', mew=0.75,
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








pltLims = ( min([dl.min() for dl in [DL1, DL2, DL3, DL4]]),
            max([dl.max() for dl in [DL1, DL2, DL3, DL4]]) )


plt.clf()
plt.figure(figsize=(8,7))
plt.axvline(zone1, 0, 1, color='green')
plt.axvline(zone2, 0, 1, color='blue')
plt.axvline(zone3, 0, 1, color='blue')
plt.axvline(zone4, 0, 1, color='blue')
plt.plot(DL4, 1*DL4+0, 'k-', alpha=0.5)
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, DL3, 'b.')
plt.plot(DL1, DL4, 'c.')
plt.xlim(*pltLims)
plt.ylim(*pltLims)
plt.show()





deltas = lambda x,y: abs(np.log10(x)-np.log10(y))


d12 = deltas(DL1, DL2)
d13 = deltas(DL1, DL3)
d14 = deltas(DL1, DL4)


plt.clf()
plt.figure(figsize=(8,7))
#plt.axvline(zone1, 0, 1, color='green')
#plt.axvline(zone2, 0, 1, color='blue')
#plt.axvline(zone3, 0, 1, color='blue')
#plt.axvline(zone4, 0, 1, color='blue')
#plt.plot(DL4, 1*DL4+0, 'k-', alpha=0.5)
plt.plot(redshifts, d12, 'r--')
plt.plot(redshifts, d13, 'b--')
plt.plot(redshifts, d14, 'c--')
#plt.xlim(*pltLims)
#plt.ylim(*pltLims)
plt.show()





#[print('%.6E  %.6E   %s'%(i,j, i==j)) for i,j in zip(DL1, DL2)]
#print(rms(abs(np.log10(DL1) - np.log10(DL2))))
#pltLims = (np.min([DL1.min(), DL2.min()]), 
#           np.min([DL1.max(), DL2.max()]))
#pltLims = (np.min([DL1.min(), DL2.min()]), 
#           np.min([DL1.max(), DL2.max()]))


plt.clf()
plt.figure(figsize=(8,7))
plt.axvline(zone1, 0, 1, color='green')
plt.axvline(zone2, 0, 1, color='blue')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.plot(DL1, DL2, 'r--')
plt.xlim(*pltLims)
plt.ylim(*pltLims)
#plt.xlim(0, 1.5E4)
#plt.ylim(0, 1.5E4)
#plt.xlim(0, zone1)
#plt.ylim(0, zone1)
plt.show()

