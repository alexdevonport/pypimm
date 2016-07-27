__author__ = 'alex'

from math import pi
import logging
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy.optimize
import warnings

from ProgBar import ProgBar
import fitutils as fit
rcParams.update({'figure.autolayout': True})
plt.style.use(['seaborn-paper'])


def fit_data(analysis):
    '''
    The univariate signal given by (timebase, signal) is fit to a
    damped sine of the form

    .. math::
        y(t) = A \exp(-g1*(t-t_0)) * cos(2 \pi f (t-t_0)) + B exp(-g_2 (t-t)0)).

    t0 is found first, and is taken as
    the time of the first zero-crossing after a user-specified amount of data to skip. A is
    found next, by shifting the last value of the signal to zero and taking the magnitude of
    the first data point. The frequency and damping are a little trickier to find, and are
    covered in their respective functions.

    Once estimates for the values are made this way, we refine the values by making a
    Lev-Mar best fit.

    :param timebase: time axis of signal
    :param signal:   voltage (y-values) of signal
    :param skip:     optionally, initial data to skip.
    :return:         dict containing amplitude, frequency, damping, and start time of the
                     damped cosine that best fits the signal.
    '''


    # First order of business, get relevant data from the analysis object
    signals = analysis.get_raw_data()
    configs = analysis.get_configs()
    setname = analysis.get_name()
    chi2lower = configs.getfloat('fit-params', 'reduced chi squared lower thresh')
    chi2upper = configs.getfloat('fit-params', 'reduced chi squared upper thresh')
    r2thresh = configs.getfloat('fit-params', 'signal fit r squared threshold')
    #sigerr = configs.getfloat('fit-params', 'estimated signal deviation')
    tstdev = configs.getfloat('fit-params', 'stdev measurement length')
    maxsnr = configs.getfloat('fit-params', 'max SNR')

    # Perform the signal fit procedure for every signal in the analysis
    # object's raw data
    r = {}  # temporary storage for fit results
    # set up the analysis progres bar
    nsigs = len(analysis.get_raw_data())
    if analysis.get_progress() is None:
        analysis.set_progress(ProgBar(msg = 'Analyzing '+setname+'...',maxn=nsigs, w=50))
        pb = analysis.get_progress()
    else:
        pb = analysis.get_progress()
        pb.restart('Analyzing '+setname+'...')

    for name, signal in signals.items():
        r[name] = {}
        timebase = analysis.get_timebase()

        # preprocess data
        timebase, signal, timebase_unf, signal_unf, signal_raw = fit.preprocess(timebase, signal, configs)
        fs = 1 / (timebase[1] - timebase[0])
        nstdev = int(fs * tstdev)
        #sigerr = np.std(signal_raw[-nstdev:])
        noise_sample = signal_raw[-nstdev:]
        sigerr = fit.noise_stdev(noise_sample)
        r[name]['noise sigma (mV)'] = sigerr
        # make amplitude estimate
        amplitude_est = np.max(signal[:25])
        # make frequency estimate
        frequency_est = fit.estimate_frequency(timebase, signal, name=name)
        r[name]['spectral peak'] = frequency_est
        # delta f / f
        r[name]['delta f / f'] = fit.spectrum_fmhw(signal, fs, name=name)
        # make damping estimate
        damping_est = fit.estimate_damping(timebase, signal, frequency_est, name=name)
        # With those estimates ready, we can try fitting the signal
        pguess = [amplitude_est, frequency_est, 1/damping_est, 0.1, 5, 0.0, 0.0]

        lbound=(-100, 0, 0, -1.5*amplitude_est, 0, -2, -2)
        ubound=(100, 3.5, 5, 1.5*amplitude_est, 5, 2, 2)
        zbounds = list(zip(lbound, ubound))

        if configs.getboolean('fit-params', 'use grid-lsq'):
            gstarts = str_to_floats(configs.get('fit-params', 'grid starts'))
            gstops = str_to_floats(configs.get('fit-params', 'grid stops'))
            glengths = str_to_floats(configs.get('fit-params', 'grid lengths'))
            gdims = list(zip(gstarts, gstops, glengths))
            #bestp = fit.grid_lsq(dsinplus_sp, timebase, signal, gdims)
        if configs.getboolean('fit-params', 'use shotgun-lsq'):
            spreads = str_to_floats(configs.get('fit-params', 'shotgun-lsq spreads'))
            bestp = pguess
            bestp = fit.shotgun_lsq(timebase, signal, dsinplus_sp,
                                       p0=pguess,
                                       spread=spreads,
                                       sigma=sigerr,
                                       maxiter=1000)
        try:
            #bestp, c2r = fit.minimize_reduced_chi2(dsinplus_sp, timebase, signal, bestp, sigma=sigerr)
            #bestp, c2r = fit.minimize_absolute(dsinplus_sp, timebase, signal, bestp, sigma=sigerr)
            #bestp, c2r = fit.minimize_lorentz(dsinplus_sp, timebase, signal, bestp, sigma=sigerr)
            worstp = [1, 1, 1, 1, 1, 1, 1]
            debounds = [(-1.5*amplitude_est, 1.5*amplitude_est),
                        (0.9*frequency_est, 1.1*frequency_est),
                        (0, 5),
                        (-1.5*amplitude_est, 1.5*amplitude_est),
                        (0, 5),
                        (-1, 1),
                        (0, 1)]
            #bestp, c2r = fit.de_lsq(dsinplus_sp, timebase, signal, debounds, sigma=sigerr)
            #bestp, c2r = fit.basin_lsq(dsinplus_sp, timebase, signal_unf, bestp, sigma=sigerr,
            #                           bounds=zbounds)
            warnings.filterwarnings('ignore')
            bestp, bestcov = scipy.optimize.curve_fit(dsinplus_sp, timebase, signal_unf,
                                                      p0=bestp)
            warnings.resetwarnings()
        except RuntimeError:
            logging.warning('FITTING FAILED FOR {}'.format(name))
            bestp = pguess
        bestfit = dsinplus_sp(timebase, *bestp)
        r2 = fit.nlcorr(dsinplus_sp, timebase, signal_unf, bestp)
        c2r = fit.redchi2(dsinplus_sp, timebase, signal_unf, bestp, sigma=sigerr)

        #confidence_limits = fit.conf_chi2(dsinplus_sp, timebase, signal, bestp, 2)
        #if r2thresh < r2 <= 1.0:
        if chi2lower <= c2r <= chi2upper and r2thresh <= r2 <= 1.0:
            r[name]['use for fit'] = True
        else:
            r[name]['use for fit'] = False
            logging.warning('Not using ' + name + ' in final analysis. chi^2: '+ str(c2r))
        #bestp = worstp
        r[name]['frequency'] = bestp[1]
        r[name]['damping']   = bestp[2]
        r[name]['amplitude'] = bestp[0]
        r[name]['interference amplitude'] = bestp[3]
        r[name]['interference damping'] = bestp[4]
        r[name]['time delay'] = bestp[5]
        r[name]['DC offset'] = bestp[6]
        r[name]['chi square'] = c2r
        r[name]['best fit SNR'] = fit.snr(signal_unf, bestfit)

        print('    {:20s}'.format(name), end=' '*2)
        for pk in bestp:
            print('{:<9.4g}'.format(pk), end=' '*2)
        print('')

        chisq = np.sum(np.multiply(signal - bestfit, signal - bestfit))
        r[name]['chi square'] = chisq
        fit.error_analysis(timebase, signal, bestfit, name=name)

        plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(timebase, signal_unf, 'b.', label='data')
        plt.plot(timebase, signal, 'r--', label='smoothed data')
        plt.plot(timebase, bestfit, 'g', label='fit')
        paramstr1 = r'r$^2$ = ' + "{0:.3f}\n".format(r2)
        paramstr2 = r'$\chi^2_\nu$ = ' + "{0:3.3f}".format(c2r)
        plt.text(0.75, 0.25,paramstr1+paramstr2,
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax.transAxes)
        plt.xlabel('time (ns)')
        plt.ylabel('signal (mV)')
        plt.title(name+' best-fit curve')
        plt.legend()
        plt.savefig(r'./sigfits/'+name+'.png')
        plt.clf()
        plt.close()
        del fig, ax
        #pb.update(1)

    # All of the fits are in r, so we'll add that to the analysis object
    analysis.set_fits(r)
    return None


def dsinplus_sp(x, p0, p1, p2, p3, p4, p5, p6):
    """
    Exactly the same as dsin, but with separate args because that's what Scipy's curve_fit takes.
    :param x:
    :param p:
    :return:
    """
    x = np.array(x)
    y =   p0 * np.cos(2*pi*p1*(x-p5)) * np.exp(-(x-p5)/abs(p2)) \
        + p3 * np.exp(-(x-p5)/abs(p4))
    return y

def str_to_floats(s):
    s = re.sub('[\[\]\(\)\,]', ' ', s)
    return [float(elt) for elt in s.split()]

def sos(x):
    return np.sum(np.multiply(x, x))
