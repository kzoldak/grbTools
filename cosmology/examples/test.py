#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 21:26:02 2020

@author: kimzoldak
"""


import numpy as np

import astropy
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import w0wzCDM
from astropy.cosmology import w0waCDM

from lumdist import Cosmology as Cosmo
from Zoldak.Math.tools import root_mean_square as rms

redshifts = np.linspace(0.1, 8, 10)

# Flat Lambda CDM: cm
cosmoconstants = {'H0':67.8, 'Om': 0.308}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='FlatLambdaCDM', 
                                 units='cm',
                                 **cosmoconstants)
cosmo = FlatLambdaCDM(H0=cosmoconstants['H0'], 
                      Om0=cosmoconstants['Om'])
#cosmo = FlatLambdaCDM(H0=65.0, 
#                      Om0=0.3)
DL2 = np.asarray(cosmo.luminosity_distance(redshifts)*u.Mpc.to(u.cm))
[print('%.6E  %.6E   %s'%(i,j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))
pltLims = (np.min([DL1.min(), DL2.min()]), 
           np.min([DL1.max(), DL2.max()]))
plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.xlim(*pltLims)
plt.ylim(*pltLims)
plt.show()


# Flat Lambda CDM: pc
cosmoconstants = {'H0':67.8, 'Om': 0.308}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='FlatLambdaCDM', 
                                 units='pc',
                                 **cosmoconstants)
cosmo = FlatLambdaCDM(H0=cosmoconstants['H0'], 
                      Om0=cosmoconstants['Om'])
#cosmo = FlatLambdaCDM(H0=70.0, 
#                      Om0=0.3)
DL2 = np.asarray(cosmo.luminosity_distance(redshifts)*u.Mpc.to(u.pc))
[print('%.6E  %.6E   %s'%(i,j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))



pltLims = (np.min([DL1.min(), DL2.min()]), 
           np.min([DL1.max(), DL2.max()]))
plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.xlim(*pltLims)
plt.ylim(*pltLims)
plt.show()




# Flat Lambda CDM: Mpc
cosmoconstants = {'H0':67.8, 'Om': 0.308}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='FlatLambdaCDM', 
                                 units='Mpc',
                                 **cosmoconstants)

cosmo = FlatLambdaCDM(H0=cosmoconstants['H0'], 
                      Om0=cosmoconstants['Om'])
#cosmo = FlatLambdaCDM(H0=70.0, 
#                      Om0=0.3)
DL2 = np.asarray(cosmo.luminosity_distance(redshifts))
[print('%.6E  %.6E   %s'%(i,j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))
pltLims = (np.min([DL1.min(), DL2.min()]), 
           np.min([DL1.max(), DL2.max()]))
plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.xlim(*pltLims)
plt.ylim(*pltLims)
plt.show()








# w0wzCDM: cm
# w0 = -1.31, wp = 1.48
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-1.31, 'wp':1.48}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0wpCDM', 
                                 units='cm',
                                 **cosmoconstants)

cosmo = w0wzCDM(H0=cosmoconstants['H0'], 
                Om0=cosmoconstants['Om'],
                Ode0=1-cosmoconstants['Om'],
                w0=cosmoconstants['w0'], 
                wz=cosmoconstants['wp'])
DL2 = np.asarray(cosmo.luminosity_distance(redshifts)*u.Mpc.to(u.cm))
[print('%.6E  %.6E   %s'%(i,j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))
pltLims = (np.min([DL1.min(), DL2.min()]), 
           np.min([DL1.max(), DL2.max()]))
plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.xlim(*pltLims)
plt.ylim(*pltLims)
plt.show()



# w0wzCDM: pc
# w0 = -1.31, wp = 1.48
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-1.31, 'wp':1.48}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0wpCDM', 
                                 units='pc',
                                 **cosmoconstants)

cosmo = w0wzCDM(H0=cosmoconstants['H0'], 
                Om0=cosmoconstants['Om'],
                Ode0=1-cosmoconstants['Om'],
                w0=cosmoconstants['w0'], 
                wz=cosmoconstants['wp'])
DL2 = np.asarray(cosmo.luminosity_distance(redshifts)*u.Mpc.to(u.pc))
[print('%.6E  %.6E   %s'%(i,j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))
pltLims = (np.min([DL1.min(), DL2.min()]), 
           np.min([DL1.max(), DL2.max()]))
plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.xlim(*pltLims)
plt.ylim(*pltLims)
plt.show()


# w0wzCDM: Mpc
# w0 = -1.31, wp = 1.48
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-1.31, 'wp':1.48}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0wpCDM', 
                                 units='Mpc',
                                 **cosmoconstants)

cosmo = w0wzCDM(H0=cosmoconstants['H0'], 
                Om0=cosmoconstants['Om'],
                Ode0=1-cosmoconstants['Om'],
                w0=cosmoconstants['w0'], 
                wz=cosmoconstants['wp'])
DL2 = np.asarray(cosmo.luminosity_distance(redshifts))
[print('%.6E  %.6E   %s'%(i, j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))
pltLims = (np.min([DL1.min(), DL2.min()]), 
           np.min([DL1.max(), DL2.max()]))
plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.xlim(*pltLims)
plt.ylim(*pltLims)
plt.show()



# w0wzCDM: Mpc
# w0 = -1.31, wp = 1.48
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-1.31, 'wp':1.48}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0wpCDM', 
                                 units='Mpc',
                                 **cosmoconstants)
cosmo = w0wzCDM(H0=cosmoconstants['H0'], 
                Om0=cosmoconstants['Om'],
                Ode0=1-cosmoconstants['Om'],
                w0=cosmoconstants['w0'], 
                wz=cosmoconstants['wp'])
DL2 = np.asarray(cosmo.luminosity_distance(redshifts))
[print('%.6E  %.6E   %s'%(i, j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))











# w0waCDM: cm
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-0.9, 'wa':0.2}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0waCDM', 
                                 units='cm',
                                 **cosmoconstants)

cosmo = w0waCDM(H0=cosmoconstants['H0'], 
                Om0=cosmoconstants['Om'],
                Ode0=1-cosmoconstants['Om'],
                w0=cosmoconstants['w0'], 
                wa=cosmoconstants['wa'])
DL2 = np.asarray(cosmo.luminosity_distance(redshifts)*u.Mpc.to(u.cm))
[print('%.6E  %.6E   %s'%(i,j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))
pltLims = (np.min([DL1.min(), DL2.min()]), 
           np.min([DL1.max(), DL2.max()]))
plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.xlim(*pltLims)
plt.ylim(*pltLims)
plt.show()






# w0waCDM: pc
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-0.9, 'wa':0.2}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0waCDM', 
                                 units='pc',
                                 **cosmoconstants)

cosmo = w0waCDM(H0=cosmoconstants['H0'], 
                Om0=cosmoconstants['Om'],
                Ode0=1-cosmoconstants['Om'],
                w0=cosmoconstants['w0'], 
                wa=cosmoconstants['wa'])
DL2 = np.asarray(cosmo.luminosity_distance(redshifts)*u.Mpc.to(u.pc))
[print('%.6E  %.6E   %s'%(i,j, i==j)) for i,j in zip(DL1, DL2)]
pltLims = (np.min([DL1.min(), DL2.min()]), 
           np.min([DL1.max(), DL2.max()]))
plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.xlim(*pltLims)
plt.ylim(*pltLims)
plt.show()



# w0waCDM: Mpc
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-0.9, 'wa':0.2}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0waCDM', 
                                 units='Mpc',
                                 **cosmoconstants)

cosmo = w0waCDM(H0=cosmoconstants['H0'], 
                Om0=cosmoconstants['Om'],
                Ode0=1-cosmoconstants['Om'],
                w0=cosmoconstants['w0'], 
                wa=cosmoconstants['wa'])
DL2 = np.asarray(cosmo.luminosity_distance(redshifts))
[print('%.6E  %.6E   %s'%(i, j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))
pltLims = (np.min([DL1.min(), DL2.min()]), 
           np.min([DL1.max(), DL2.max()]))
plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.xlim(*pltLims)
plt.ylim(*pltLims)
plt.show()



import matplotlib.pyplot as plt
# %matplotlib inline


redshifts = np.linspace(0, 10, 100)


# w0waCDM vs w0wzCDM: Mpc
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-0.9, 'wa':0.2}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0waCDM', 
                                 units='Mpc',
                                 **cosmoconstants)
zone1 = cosmo.luminosity_distance(redshifts=[1.0,], 
                                 model='w0waCDM', 
                                 units='Mpc',
                                 **cosmoconstants)
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-1.31, 'wp':1.48}
DL2 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0wpCDM', 
                                 units='Mpc',
                                 **cosmoconstants)
zone2 = cosmo.luminosity_distance(redshifts=[1.0,], 
                                 model='w0wpCDM', 
                                 units='Mpc',
                                 **cosmoconstants)
#[print('%.6E  %.6E   %s'%(i, j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))
pltLims = (np.min([DL1.min(), DL2.min()]), 
           np.min([DL1.max(), DL2.max()]))

plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r-')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.axvline(zone1, 0, 1, color='green')
plt.axvline(zone2, 0, 1, color='blue')
#plt.xlim(*pltLims)
#plt.ylim(*pltLims)
#plt.xlim(0, 1.5E4)
#plt.ylim(0, 1.5E4)
plt.xlim(0, zone1)
plt.ylim(0, zone1)
plt.show()






# w0waCDM vs w0wzCDM: Mpc
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-0.9, 'wa':0.2}
cosmo = Cosmo()
DL1 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0waCDM', 
                                 units='cm',
                                 **cosmoconstants)
zone1 = cosmo.luminosity_distance(redshifts=[1.0,], 
                                 model='w0waCDM', 
                                 units='cm',
                                 **cosmoconstants)
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-1.31, 'wp':1.48}
DL2 = cosmo.luminosity_distance(redshifts=redshifts, 
                                 model='w0wpCDM', 
                                 units='cm',
                                 **cosmoconstants)
zone2 = cosmo.luminosity_distance(redshifts=[1.0,], 
                                 model='w0wpCDM', 
                                 units='cm',
                                 **cosmoconstants)
[print('%.6E  %.6E   %s'%(i, j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))
pltLims = (np.min([DL1.min(), DL2.min()]), 
           np.min([DL1.max(), DL2.max()]))
plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.axvline(zone1, 0, 1, color='green')
plt.axvline(zone2, 0, 1, color='blue')
plt.xlim(*pltLims)
plt.ylim(*pltLims)
plt.show()






cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-0.9, 'wa':0.2}
cosmo = w0waCDM(H0=cosmoconstants['H0'], 
                Om0=cosmoconstants['Om'],
                Ode0=1-cosmoconstants['Om'],
                w0=cosmoconstants['w0'], 
                wa=cosmoconstants['wa'])
cosmoconstants = {'H0':67.8, 'Om': 0.308, 'w0':-1.31, 'wp':1.48}
cosmo = w0wzCDM(H0=cosmoconstants['H0'], 
                Om0=cosmoconstants['Om'],
                Ode0=1-cosmoconstants['Om'],
                w0=cosmoconstants['w0'], 
                wz=cosmoconstants['wp'])
DL2 = np.asarray(cosmo.luminosity_distance(redshifts))
[print('%.6E  %.6E   %s'%(i, j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))

plt.clf()
plt.figure(figsize=(8,7))
plt.plot(DL1, DL2, 'r.')
plt.plot(DL1, 1*DL1+0, 'k-', alpha=0.5)
plt.xlim(0, 1.1E5)
plt.ylim(0, 1.1E5)
plt.show()



DL2 = np.asarray(cosmo.luminosity_distance(redshifts))
[print('%.6E  %.6E   %s'%(i, j, i==j)) for i,j in zip(DL1, DL2)]
print(rms(abs(np.log10(DL1) - np.log10(DL2))))
