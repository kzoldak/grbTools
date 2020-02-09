#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 20:42:24 2020

@author: kimzoldak
"""

import numpy as np


import matplotlib.pyplot as plt
%matplotlib inline


# Choose the "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534


# Generate some synthetic data from the model.
N = 50
x = np.sort(10*np.random.rand(N))
yerr = 0.1+0.5*np.random.rand(N)
y = m_true*x+b_true
y += np.abs(f_true*y) * np.random.randn(N)
y += yerr * np.random.randn(N)


A = np.vstack((np.ones_like(x), x)).T
C = np.diag(yerr * yerr)
cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))



def lnlike(theta, x, y, yerr):
    m, b, lnf = theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))



import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
m_ml, b_ml, lnf_ml = result["x"]



def lnprior(theta):
    m, b, lnf = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
        return 0.0
    return -np.inf


def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)


ndim, nwalkers = 3, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]


import emcee

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))

sampler.run_mcmc(pos, 500)


# lengths: [100][500][3]
#  sampler.chain[99][499][2]

pltKwgs = dict(color='k', alpha=0.1)
labels = ['m', 'b', 'f']
plt.clf()
for i in range(sampler.chain.shape[2]):
    plt.figure(figsize=(11,5))
    [plt.plot(sampler.chain[j][:,0], **pltKwgs) for j in range(0, 100)]
    plt.ylabel(labels[i])
    plt.xlabel('step number')
    plt.xlim(0, 500)
    plt.show()

samples = sampler.chain[:, 50:, :].reshape((-1, ndim))



pltKwgs = dict(color='k', alpha=0.25)
labels = ['m', 'b', 'f']
plt.clf()
for i in range(samples.shape[1]):
    plt.figure(figsize=(11,5))
    plt.plot(samples[:,0], **pltKwgs)
    plt.ylabel(labels[i])
    plt.xlabel('step number')
    plt.xlim(0, samples.shape[0])
    plt.show()


import corner
fig = corner.corner(samples, labels=["$m$", "$b$", "$\ln\,f$"],
                      truths=[m_true, b_true, np.log(f_true)])


#fig.savefig("triangle.png")


xl = np.array([0, 10])
for m,b,lnf in samples[np.random.randint(len(samples), size=100)]:
    plt.plot(xl, m*xl+b, color="k", alpha=0.1)
plt.plot(xl, m_true*xl+b_true, color="r", lw=2, alpha=0.8)
plt.errorbar(x, y, yerr=yerr, fmt=".k")




samples[:, 2] = np.exp(samples[:, 2])
m_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))


print(m_mcmc)
print(b_mcmc)
print(f_mcmc)
