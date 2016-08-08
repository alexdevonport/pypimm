import scipy.signal
import numpy as np
import numpy.fft as ft
import scipy.optimize
import matplotlib.pyplot as plt
import os
import fitutils as fit
import random

__author__ = 'alex'

def optimal_filter(y, fs, fmodel, f0, nmodel, n0, name=''):

    mf, mn = np.size(f0), np.size(n0)
    yf = ft.rfft(y)
    f = fs * np.arange(np.size(yf)) / np.size(y)

    def tmodel(f1, *p0):
        p1 = p0
        pn = p1[:mn]
        pf = p1[mn:]
        return np.add(fmodel(f1, *pf), nmodel(f1, *pn))
    pguess = f0 + n0
    bestp, _ = fit.basin_lsq(tmodel, f, np.abs(yf), pguess, niter=25)
    bestp, _ = fit.minimize_lorentz(tmodel, f, np.abs(yf), pguess)
    #bestp, pcov = scipy.optimize.curve_fit(tmodel, f, np.abs(yf), pguess)
    r2 = fit.nlcorr(tmodel, f, np.abs(yf), bestp)
    #print('OPTIMAL FILTER R SQUARE: ', r2)
    bestfit = tmodel(f, *bestp)
    bestpn = bestp[:mn]
    bestpf = bestp[mn:]
    #s = np.power(np.abs(fmodel(f, *bestpf)), 2)
    n = np.power(np.abs(nmodel(f, *bestpn)), 2)
    c = np.power(np.abs(yf), 2)
    smodel = np.abs(yf) - nmodel(f, *bestpn)
    s = np.power(np.abs(smodel), 2)
    #s = c - n
    phi = np.divide(s, np.add(s, n))
    yfphi = np.multiply(yf, phi)
    yphi = ft.irfft(yfphi)

    pltname = str(random.randint(1, 100000000))
    plt.clf()
    plt.plot(f, np.abs(yf), label='data')
    plt.plot(f, bestfit, label='fit')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('PSD (V$^2$ / GHz)')
    plt.title(name + ' Spectrum')
    fp = os.path.join('.','spectra', pltname + ' spectrum.png')
    plt.savefig(fp)

    plt.clf()
    plt.plot(yphi)
    plt.title(name + ' Spectrum')
    fp = os.path.join('.','spectra', pltname + ' filtered.png')
    plt.savefig(fp)

    if r2 < 0.1:
        raise RuntimeError('Could not construct optimal filter')

    return yphi