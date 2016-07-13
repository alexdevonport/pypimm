__author__ = 'alex'

from math import pi
import logging

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

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
    r2thresh = analysis.get_configs('signal fit r squared threshold')

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
        timebase, signal = fit.preprocess(timebase, signal, configs)
        fs = 1 / (timebase[1] - timebase[0])
        # make amplitude estimate
        amplitude_est = np.max(signal[:25])
        # make frequency estimate
        frequency_est = fit.estimate_frequency(timebase, signal, name=name)
        # delta f / f
        r[name]['delta f / f'] = fit.spectrum_fmhw(signal, fs, name=name)
        # make damping estimate
        damping_est = fit.estimate_damping(timebase, signal, frequency_est, name=name)
        # With those estimates ready, we can try fitting the signal
        try:
            sigpopt, sigpcov, r2 = fit.shotgun_lsq(timebase, signal, dsinplus_sp,
                                           p0=[amplitude_est, frequency_est, 1/damping_est, 0.3, 1, 0.1],
                                           spread=[1.0, 0.4, 3, 0.4, 2, 0.1],
                                           sigma=0.1,
                                           maxiter=1000)
            if r2 < r2thresh:
                r[name]['use for fit'] = False
                logging.warning('Not using ' + name + ' in final analysis. r^2: '+ str(r2))
            else:
                r[name]['use for fit'] = True
            r[name]['frequency'] = sigpopt[1]
            r[name]['damping']   = sigpopt[2]
            r[name]['amplitude'] = sigpopt[0]
            r[name]['chi square'] = r
            bestfit = dsinplus_sp(timebase, *sigpopt)
            chisq = np.sum(np.multiply(signal - bestfit, signal - bestfit))
            r[name]['chi square'] = chisq
        except RuntimeError:
            r2 = 0
            r[name]['use for fit'] = False
            r['chi square'] = 0
            #print('Fitting failed for signal '+name+'. Using estimate values.')
            r[name]['frequency'] = frequency_est
            r[name]['damping']   = damping_est
            r[name]['amplitude'] = amplitude_est
            bestfit = dsinplus_sp(timebase, *[amplitude_est, frequency_est, 1/damping_est, 0.3, 1, 0.1])

        fit.error_analysis(timebase, signal, bestfit, name=name)

        plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(timebase, signal, 'b.', label='data')
        plt.plot(timebase, bestfit, 'g', label='fit')
        paramstr = 'r$^2$ = ' + "{0:.3f}".format(r2)
        plt.text(0.75, 0.25,paramstr,
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax.transAxes)
        plt.xlabel('time (ns)')
        plt.ylabel('signal (V)')
        plt.title(name+' best-fit curve')
        plt.legend()
        plt.savefig(r'./sigfits/'+name+'.png')
        plt.clf()
        plt.close()
        del fig, ax
        pb.update(1)

    # All of the fits are in r, so we'll add that to the analysis object
    analysis.set_fits(r)
    return None


def dsinplus_sp(x, p0, p1, p2, p3, p4, p5):
    """
    Exactly the same as dsin, but with separate args because that's what Scipy's curve_fit takes.
    :param x:
    :param p:
    :return:
    """
    x = np.array(x)
    y =   p0 * np.cos((x-p5)*2*pi*p1)*np.exp(-(x-p5)/abs(p2)) \
        + p3 * np.exp(-(x-p5)/abs(p4))
    y[y == np.inf] = 0
    y[y == np.nan] = 0
    return y
