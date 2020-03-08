#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 13:35:12 2020

@author: kimzoldak
"""

import os
from astropy.io import fits as pyfits

datadir = '/Users/kimzoldak/RMFITResults/bn080916009/integrated/'
fn = 'rmfit_080916C_int_band_GBM.fit'

f = open(os.path.join(datadir, fn))


class RMFIT(object):
    def __init__(self, modelnames):
        self.modelnames = ['band', 'sbpl', 'copl', 'lpow', 'bbody', ]