__author__ = 'alex'

from math import pi
import logging
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy.optimize
import warnings
import collections
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
    fields = analysis.get_fields()
    setname = analysis.get_name()
    chi2lower = configs.getfloat('fit-params', 'reduced chi squared lower thresh')
    chi2upper = configs.getfloat('fit-params', 'reduced chi squared upper thresh')
    r2thresh = configs.getfloat('fit-params', 'signal fit r squared threshold')
    tstdev = configs.getfloat('fit-params', 'stdev measurement length')

    # Perform the signal fit procedure for every signal in the analysis
    # object's raw data
    r = collections.OrderedDict()  # temporary storage for fit results
    # set up the analysis progres bar
    nsigs = len(analysis.get_raw_data())
    if analysis.get_progress() is None:
        analysis.set_progress(ProgBar(msg = 'Analyzing '+setname+'...',maxn=nsigs, w=50))
        pb = analysis.get_progress()
    else:
        pb = analysis.get_progress()
        pb.restart('Analyzing '+setname+'...')

    # print header for tabulated printout
    print('    {:20s}'.format(' '), end=' '*2)
    for pk in ['ampl (mV)', 'freq(GHz)', 'tau', 'intampl', 'inttau', 'del (ns)', 'offs (mV)']:
            print('{:9s}'.format(pk), end=' '*2)
    print('')
    for (name, signal), h in zip(signals.items(), fields):
        fitres = collections.OrderedDict()
        fitres['bias field'] = h
        timebase = analysis.get_timebase()
        timebase, signal, timebase_unf, signal_unf, signal_raw = fit.preprocess(timebase, signal, configs)
        fs = 1 / (timebase[1] - timebase[0])
        nstdev = int(fs * tstdev)
        noise_sample = signal_raw[-nstdev:]
        fnoise_sample = signal[-nstdev:]
        sigerr = fit.noise_stdev(noise_sample)
        fsigerr = fit.noise_stdev(fnoise_sample)
        fitres['noise sigma (mV)'] = sigerr
        amplitude_est = np.max(signal[:25])
        frequency_est = fit.estimate_frequency(timebase, signal_unf, name=name)
        fitres['spectral peak'] = frequency_est
        fitres['delta f / f'] = fit.spectrum_fmhw(signal, fs, name=name)
        damping_est = fit.estimate_damping(timebase, signal, frequency_est, name=name)
        # With those estimates ready, we can try fitting the signal
        pguess = [amplitude_est, frequency_est, 1/damping_est, 0.1, 5, 0.0, 0.0]

        lbound=(-100, 0, 0, -1.5*amplitude_est, 0, -2, -100)
        ubound=(100, 3.5, 5, 1.5*amplitude_est, 5, 2, 100)
        zbounds = list(zip(ubound, lbound))

        if configs.getboolean('fit-params', 'use shotgun-lsq'):
            spreads = str_to_floats(configs.get('fit-params', 'shotgun-lsq spreads'))
            bestp = pguess
            bestp = fit.shotgun_lsq(timebase, signal, sfit,
                                       p0=pguess,
                                       spread=spreads,
                                       sigma=fsigerr,
                                       maxiter=1000)
        try:
            bestp, c2r = fit.basin_lsq(sfit, timebase, signal, bestp, sigma=fsigerr,
                                       bounds=zbounds)
            warnings.filterwarnings('ignore')
            bestp, bestcov = scipy.optimize.curve_fit(sfit, timebase, signal_unf,
                                                      sigma=sigerr,
                                                      absolute_sigma=True,
                                                      method='trf',
                                                      p0=bestp)
            warnings.resetwarnings()
            #print('BEFORE', bestp)
            #worstint = fit.conf_chi2(sfit, timebase, signal_unf, bestp, 2, sigma=sigerr)
            #print('AFTER', bestp)
            bestpc = fit.conf1d(bestp, bestcov, 2)
            #bestint = [st[1] for st in bestpc]
            #print('DELTA CONFIDENCE INTERVALS: ', np.subtract(bestint, worstint))

        except RuntimeError:
            logging.warning('FITTING FAILED FOR {}'.format(name))
            bestp = pguess
            bestpc =fit.conf1d(bestp, np.zeros((len(bestp), len(bestp))))
        bestfit = sfit(timebase, *bestp)
        r2 = fit.nlcorr(sfit, timebase, signal_unf, bestp)
        c2r = fit.redchi2(sfit, timebase, signal_unf, bestp, sigma=sigerr)
        fit.error_analysis(timebase, signal_unf, bestfit, name=name)
        #if r2thresh < r2 <= 1.0:
        if chi2lower <= c2r <= chi2upper and r2thresh <= r2 <= 1.0:
            use_for_fit = True
        else:
            logging.warning('Not using ' + name + ' in final analysis. chi^2: '+ str(c2r))
            use_for_fit = False

        # give names to fit parameters, for convenience in excel spreadsheet
        fitres['frequency'] = np.abs(bestpc[1][0])
        fitres['frequency interval'] = bestpc[1][1]
        fitres['frequency sigma'] = bestpc[1][2]
        fitres['damping']   = bestpc[2][0]
        fitres['damping interval']   = bestpc[2][1]
        fitres['amplitude'] = bestpc[0][0]
        fitres['amplitude interval'] = bestpc[0][1]
        fitres['interference amplitude'] = bestpc[3][0]
        fitres['interference damping'] = bestpc[4][0]
        fitres['time delay'] = bestpc[5][0]
        fitres['DC offset'] = bestpc[6][0]
        fitres['chi square'] = c2r
        fitres['best fit SNR'] = fit.snr(signal_unf, bestfit)

        # If this signal matches the goodness-of-fit criteria,
        # add it to the fit results dictionary
        if use_for_fit:
            r[name] = fitres

        print('    {:20s}'.format(name), end=' '*2)
        for pk in bestp:
            print('{:<9.4g}'.format(pk), end=' '*2)
        print('')

        plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(timebase, signal_unf, 'b.', label='data')
        plt.plot(timebase, signal, 'r--', label='smoothed data')
        plt.plot(timebase, bestfit, 'g', label='fit')
        paramstr1 = r'r$^2$ = ' + "{0:.3f}\n".format(r2)
        #paramstr1 = r''
        paramstr2 = r'$\chi^2_\nu$ = ' + "{0:3.3f}\n".format(c2r)
        paramstr3 = r'$F_p$ = {:.3g} $\pm$ {:.2g} GHz'.format(bestpc[1][0],bestpc[1][1])
        plt.text(0.75, 0.25,paramstr1+paramstr2+paramstr3,
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


def sfit(x, p0, p1, p2, p3, p4, p5, p6):
    x = np.array(x)
    y =   (p0 * np.cos(2*pi*p1*(x-p5)) * np.exp(-(x-p5)*abs(p2))
        + p3 * np.exp(-(x-p5)/abs(p4)) + p6)
    return y

def msfit(x, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9):
    x = np.array(x)
    y = (p0 * np.cos(2*pi*p1*(x-p2)) * np.exp(-(x-p2)*abs(p3))
        + p4 * np.cos(2*pi*p5*(x-p2)) * np.exp(-(x-p2)*abs(p6))
        + p7 * np.exp(-(x-p2)/abs(p8)) + p9)
    return y


def str_to_floats(s):
    s = re.sub('[\[\]\(\)\,]', ' ', s)
    return [float(elt) for elt in s.split()]

def sos(x):
    return np.sum(np.multiply(x, x))
