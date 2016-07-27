import re
from math import pi
import os
import scipy.optimize
import matplotlib.pyplot as plt
import numpy as np
import fitutils as fit
from collections import OrderedDict

__author__ = 'alex'


# def pimm_characterize_fits(fits, name=None):
def characterize_fits(analysis):
    '''

    :param fits:
    :return:
    '''

    # get relevant data from the analysis object
    fits = analysis.get_fits()
    name = analysis.get_name()
    configs = analysis.get_configs()

    # placeholder for characterisation results
    r = OrderedDict()

    # get an ordered list of all bias fields
    # TODO: remove the neeed for this by using ordered dicts
    key_pattern = re.compile('([+-]?\d+\.?\d*)(.*)')
    hs = []
    amps = []
    fs = []
    ds = []
    cs = []
    dffs = []
    sps = []
    snrs = []

    # Some of the calculated values in fs and ds may be bad, if we couldn't fit the
    # signal. If that's the case, those signals will have been flagged. We'll get
    # rid of them now.

    for key in fits.keys():
        mo = key_pattern.match(key)
        if not mo == None:
            if fits[key]['use for fit']:
                h = float(mo.group(1))
                hs.append(h)
                fs.append(fits[key]['frequency'])
                ds.append(fits[key]['damping'])
                amps.append(fits[key]['amplitude'])
                cs.append(fits[key]['chi square'])
                dffs.append(fits[key]['delta f / f'])
                sps.append(fits[key]['spectral peak'])
                snrs.append(fits[key]['best fit SNR'])

    # Since the dict was unordered this whole time, we need to sort all
    # three of these together
    hs, fs, ds, amps, cs, dffs, sps, snrs = (list(t) for t in zip(*sorted(zip(hs, fs, ds, amps, cs, dffs, sps, snrs))))
    r['Bias field (Oe)'] = hs
    r['Precessional frequency (GHz)'] = fs
    r['amplitude'] = amps
    r['Chi square'] = cs
    r['delta f / f'] = dffs
    r['spectral peak'] = sps
    r['SNR'] = snrs

    # The damping values may be negative. This is because the absolute value is
    # used in the fit to remove negative guesses. We take care of that here.
    # Very rarely, the same thing happens for frequency, too.
    ds = np.abs(ds)
    fs = np.abs(fs)

    # Prepare H and f data for calculating Ms and Hk
    # TODO: remove all frequency data points below a certain field, so that Hk gets
    # TODO: fit correctly
    hkthresh = configs.getfloat('characterize', 'hk field fit thresh')
    hs_msfit = []
    fs_msfit = []
    sps_msfit = []
    for h, f, s in zip(hs, fs, sps):
        if abs(h) > hkthresh:
            hs_msfit.append(h)
            fs_msfit.append(f)
            sps_msfit.append(s)

    hsi = 1000 / (4 * pi) * np.array(hs)  # H, in SI units (A/M)
    hsi_msfit = 1000 / (4 * pi) * np.array(hs_msfit)  # H, in SI units (A/M)
    omegas = 2 * pi * np.array(fs_msfit) * 1E9
    omegas_sps = 2 * pi * np.array(sps) * 1E9

    hkguess = 10 * 1000 / (4 * pi)
    msguess = 800 * 1E3
    hcpguess = 1 * 1000 / (4 * pi)

    #uncomment these lines to use spectral peaks for mhsk estimate
    #omegas = np.array(sps_msfit) * 2 * pi * 1E9

#    fitbounds = [
#        (0,2000*1E3),
#        (0, 75 * 1000 / (4 * pi)),
#        (-10 * 1000 / (4 * pi), 10 * 1000 / (4 * pi))
#    ]
    fitbounds = ([0, 0, -10 * 1000/(4*pi)],
                 [2000*1E3, 75*1000/(4*pi), 10*1000/(4*pi)])

    # fpopt, fpcov = sp.optimize.curve_fit(precession, hsi_msfit, omegas, p0=[msguess, hkguess, hcpguess])
    #bestp = fit.shotgun_lsq(hsi_msfit, omegas, precession,
    #                        p0=[msguess, hkguess, hcpguess],
    #                        spread=[100 * 1E3, 2 * 1000 / (4 * pi), 0.5 * 1000 / (4 * pi)], sigma=1, maxiter=10000)

    #try:
    bestp, bestcov = scipy.optimize.curve_fit(precession, hsi_msfit,
                                              omegas, p0=[msguess, hkguess, hcpguess],
                                              bounds=fitbounds)
#except:
    #    print('Could not characterize Ms and Hk.')
    #    bestp, bestcov = [1,1,1], []
    bestp, _ = fit.minimize_lorentz(precession, hsi_msfit, omegas, bestp, sigma=0.1E9)
    print(bestp)
    bestfit = precession(hsi_msfit, *bestp)
    ssres = sos(omegas - bestfit)  # sum of squares of the residuals
    sigmean = np.mean(omegas)
    sstot = sos(omegas - sigmean)  # total sum of squares
    r2 = 1.0 - ssres / sstot

    hs_bestfit = np.linspace(np.min(hs), np.max(hs), 1000)
    hsi_bestfit = hs_bestfit * 1000 / (4 * pi)
    omega_bestfit = precession(hsi_bestfit, *bestp)
    mscgs = bestp[0] * 1E-3
    hkcgs = bestp[1] * 4 * pi / 1000
    hcpcgs = bestp[2] * 4 * pi / 1000

    # prepare frequency vs field plot
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(hs_msfit, omegas, 'bo ', label='data')
    plt.plot(hs_bestfit, omega_bestfit, label='fit')
    plt.legend()
    plt.xlabel('bias field (Oe)')
    plt.ylabel('$\omega_p$ (rad/s)')
    plt.title('Precessional frequency\n' + name)
    paramstr = ('$M_s$ = ' + "{0:.2f}".format(mscgs) + ' emu/cm^3\n'
                + '$H_k$ = ' + "{0:.2f}".format(hkcgs) + ' Oe\n'
                + '$H_{cp}$ = ' + "{0:.2f}".format(hcpcgs) + ' Oe\n'
                + '$r^2$ = ' + "{0:.4f}".format(r2))
    plt.text(0.75, 0.25, paramstr,
             horizontalalignment='center',
             verticalalignment='center',
             transform=ax.transAxes)
    fp = os.path.join('.', name + '-f-vs-h.png')
    plt.savefig(fp, dpi=150)
    plt.clf()
    r['Ms (emu/cm^3)'] = mscgs
    r['Hk (Oe)'] = hkcgs

    # TODO: make a function to clean up plotting section a little

    # prepare damping vs field plot
    gmr = 28 * 2 * pi  # Gyromagnetic ratio, GHz/T
    mu0 = 4 * pi * 1E-7  # vacuum permeability (T-m/A)
    ds = np.divide(2, ds) * 1 / (gmr * mu0 * mscgs * 1E3)
    r['average damping'] = np.mean(ds)
    r['Damping'] = ds

    plt.plot(hs, ds, 'bs-')
    plt.xlabel('bias field (Oe)')
    plt.ylabel('damping')
    plt.title('Damping\n' + name[:10] + '...')
    plt.savefig(r'./' + name + '-d-vs-h.png')
    plt.clf()

    # prepare amplitude vs field plot
    plt.plot(hs, amps, 'bo ')
    plt.xlabel('bias field (Oe)')
    plt.ylabel('amplitude (mV)')
    plt.title('sinusoidal fit amplitude\n' + name[:10] + '...')
    plt.savefig(os.path.join('.', name + '-ampl.png'))
    plt.clf()


    # prepare chi square vs H plot and chi square histogram
    plt.plot(hs, cs)
    plt.xlabel('bias field (Oe)')
    plt.ylabel('Chi square')
    plt.title('Chi-square error of signal fits')
    plt.savefig(os.path.join('.', name + '-chi-square.png'))
    plt.clf()
    n, bins, patches = plt.hist(cs, 100, normed=1, alpha=0.75)
    plt.title('Chi-square error of signal fits')
    plt.savefig(os.path.join('.', name + '-chi-square-hist.png'))
    plt.clf()

    # Now that that's all done, add the results to the analysis object
    analysis.set_results(r)
    return None


def precession(x, p0, p1, p2):
    mub = 9.274E-24  # Bohr Magneton (J/T)
    mu0 = 4 * pi * 1E-7  # vacuum permeability (T-m/A)
    hbar = 1.06E-34  # reduced Planck's constant (J-S)
    g = 2  # Spectroscopic splitting factor, from Silva et al
    gamma = g * mub / hbar
    # fp = gamma * mu0 * np.sqrt(p0*np.abs(x + p1))
    fp = gamma * mu0 * np.sqrt(p0 * (p1 + np.abs(x + p2)))
    return fp


def sos(x):
    return np.sum(np.multiply(x, x))
