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

    #

    # placeholder for characterisation results
    r = OrderedDict()

    # get an ordered list of all bias fields
    # TODO: remove the neeed for this by using ordered dicts
    key_pattern = re.compile('([+-]?\d+\.?\d*)(.*)')
    hs = get_subkey(fits, 'bias field')
    amps, ampcs = get_subkey(fits, 'amplitude'), get_subkey(fits, 'amplitude interval')
    fs, fcs = get_subkey(fits, 'frequency'), get_subkey(fits, 'frequency interval')
    fsig = get_subkey(fits, 'frequency sigma')
    ds, dcs = get_subkey(fits, 'damping'), get_subkey(fits, 'damping interval')

    """
    for key in fits.keys():
        mo = key_pattern.match(key)
        if not mo == None:
            if fits[key]['use for fit']:
                h = float(mo.group(1))
                hs.append(h)
                fs.append(fits[key]['frequency'][0])
                fcs.append(fits[key]['frequency'][1])
                ds.append(fits[key]['damping'][0])
                dcs.append(fits[key]['damping'][1])
                amps.append(fits[key]['amplitude'][0])
                ampcs.append(fits[key]['amplitude'][1])
    """

    # The damping values may be negative. This is because the absolute value is
    # used in the fit to remove negative guesses. We take care of that here.
    # Very rarely, the same thing happens for frequency, too.
    ds = np.abs(ds)
    fs = np.abs(fs)

    # Prepare H and f data for calculating Ms and Hk
    hkthresh = configs.getfloat('characterize', 'hk field fit thresh')
    hs_msfit = []
    fs_msfit = []
    fsig_msfit = []
    fcs_msfit = []
    for h, f, s, c in zip(hs, fs, fsig, fcs):
        if abs(h) >= hkthresh:
            hs_msfit.append(h)
            fs_msfit.append(f)
            fsig_msfit.append(s)
            fcs_msfit.append(c)

    hsi_msfit = 1000 / (4 * pi) * np.array(hs_msfit)  # H, in SI units (A/M)
    omegas = 2 * pi * np.array(fs_msfit)
    omega_sig = 2 * pi * np.array(fsig_msfit)
    omegacs = 0.5 * 2 * pi * np.array(fcs_msfit)
    hkguess = 10 * 1000 / (4 * pi)
    msguess = 800 * 1E3
    hcpguess = 1 * 1000 / (4 * pi)
    r['hs'] = hs
    fitbounds = ([0, 0, -10 * 1000/(4*pi)],
                 [2000*1E3, 75*1000/(4*pi), 10*1000/(4*pi)])

    omega_sigz = []
    omegacz = []
    sigavg = np.median(omega_sig)
    cavg = np.median(omegacs)
    for sig, c in zip(omega_sig, omegacs):
        if sig < 1E-8 or sig > 1E3:
            print('conf1d FOUND A BAD CONFINT')
            omega_sigz.append(sigavg)
            omegacz.append(cavg)
        else:
            omega_sigz.append(sig)
            omegacz.append(c)

    bestp, bestcov = scipy.optimize.curve_fit(precession, hsi_msfit,
                                              omegas, p0=[msguess, hkguess, hcpguess],
                                              bounds=fitbounds, sigma=np.array(omega_sigz))

    c2r = fit.redchi2(precession, hsi_msfit, omegas,
                      bestp, sigma=omega_sig)

    bestpc = fit.conf1d(bestp, bestcov, stdevs=2)
    msconfint = 2 * bestpc[0][1] / 1000
    hkconfint = 2 * bestpc[1][1] * 4 * pi / 1000
    hcconfint = 2 * bestpc[2][1] * 4 * pi / 1000

    bestfit = precession(hsi_msfit, *bestp)
    fit.error_analysis(hsi_msfit, hsi_msfit, bestfit, name='properties')

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
    plt.plot(hs_msfit, omegas/ (2*pi), 'bo ', label='data')

    # TODO: get error bars working. Holy cow, they do not want to work.
    plt.errorbar(hs_msfit, omegas/(2*pi), yerr=np.divide(omegacz, 2*pi), fmt='o')
    plt.plot(hs_bestfit, omega_bestfit/(2*pi), label='fit')
    plt.legend()
    plt.xlabel('bias field (Oe)')
    plt.ylabel('$f_p$ (GHz)')
    plt.title('Precessional frequency\n' + name)
    paramstr = ('$M_s$ = ' + "{0:.3g}".format(mscgs)
                + ' $\pm$ {0:.3g}'.format(msconfint) + ' emu/cm^3\n'
                + '$H_k$ = ' + "{0:.2g}".format(hkcgs)
                + ' $\pm$ {0:.2g}'.format(hkconfint) + ' Oe\n'
                + '$H_{cp}$ = ' + "{0:.2f}".format(hcpcgs)
                + ' $\pm$ {0:.2g}'.format(hcconfint) + ' Oe\n'
#                + '$r^2$ = ' + "{0:.4f}".format(r2) + '\n'
                + r'$\chi^2_\nu$ = {0:.4f}'.format(c2r))
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
    ds = np.multiply(2, ds) * 1 / (gmr * mu0 * mscgs * 1E3)
    dcs = np.multiply(2, dcs) * 1 / (gmr * mu0 * mscgs * 1E3)
    r['average damping'] = np.mean(ds)
    r['Damping'] = ds

    deb = np.multiply(0.5, dcs)
    plt.plot(hs, ds, 'bs-')
    #plt.errorbar(hs, ds, yerr= deb)
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

    # Now that that's all done, add the results to the analysis object
    analysis.set_results(r)
    return None


def precession(x, ms, hk, hcp):
    mub = 9.274E-24  # Bohr Magneton (J/T)
    mu0 = 4 * pi * 1E-7  # vacuum permeability (T-m/A)
    hbar = 1.06E-34  # reduced Planck's constant (J-S)
    g = 2  # Spectroscopic splitting factor, from Silva et al
    gamma = 2 * pi * 28
    # fp = gamma * mu0 * np.sqrt(p0*np.abs(x + p1))
    fp = gamma * mu0 * np.sqrt(ms * (hk + np.abs(x + hcp)))
    return fp

def precession2(h, ms, hk, hcp):
    mu0 = 4 * pi * 1E-7  # vacuum  (T-m/A)
    gamma = 2 * pi * 28 # GHz / T
    heff = np.abs(h + hcp)
    fp = gamma * mu0 * np.sqrt((heff + hk)*(heff + hk + ms))
    return fp


def sos(x):
    return np.sum(np.multiply(x, x))

def get_subkey(d, subkey):
    r = []
    for key, val in d.items():
        try:
            r.append(val[subkey])
        except (TypeError, KeyError):
            pass
    return r

def replinf(x, repl):
    """
    removed infinite values from array x and replace them with repl.
    :param x:
    :return:
    """
    for i,_ in enumerate(x):
        if x[i] == np.inf:
            x[i] = repl
    return x