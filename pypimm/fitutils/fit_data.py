__author__ = 'alex'

from math import pi
import logging
import re
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy.optimize
import scipy.signal
import warnings
import collections
from ProgBar import ProgBar
import fitutils as fit

rcParams.update({'figure.autolayout': True})
plt.style.use(['seaborn-paper'])


def fit_data(analysis):

    # First order of business, get relevant data from the analysis object
    signals = analysis.get_raw_data()
    configs = analysis.get_configs()
    fields = analysis.get_fields()
    setname = analysis.get_name()
    chi2lower = configs.getfloat('fit-params', 
        'reduced chi squared lower thresh')
    chi2uppers = configs.getfloat('fit-params',
        'reduced chi squared upper thresh')
    chi2uppere = configs.getfloat('fit-params',
        'reduced chi squared upper thresh for expfit')
    r2thresh = configs.getfloat('fit-params',
        'signal fit r squared threshold')
    tstdev = configs.getfloat('fit-params',
        'stdev measurement length')

    # Perform the signal fit procedure for every signal in the analysis
    # object's raw data
    r = collections.OrderedDict()  # temporary storage for fit results
    rall = collections.OrderedDict()
    # set up the analysis progres bar
    nsigs = len(analysis.get_raw_data())
    if analysis.get_progress() is None:
        analysis.set_progress(ProgBar(msg = 'Analyzing '+setname+'...',
            maxn=nsigs, w=50))
        pb = analysis.get_progress()
    else:
        pb = analysis.get_progress()
        pb.restart('Analyzing '+setname+'...')

    # print header for tabulated printout
    print('    {:20s}'.format(' '), end=' '*2)
    for pk in ['ampl (mV)', 'freq(GHz)', 'tau', 'chi^2']:
            print('{:9s}'.format(pk), end=' '*2)
    print('')
    for (name, signal), h in zip(signals.items(), fields):
        fitres = collections.OrderedDict()
        fitres['bias field'] = h
        timebase = analysis.get_timebase()
        timebase, signal, timebase_unf, signal_unf, signal_raw, tpeak = fit.preprocess(timebase, signal, configs)
        fs = 1 / (timebase[1] - timebase[0])
        nstdev = int(fs * tstdev)
        noise_sample = signal_raw[-nstdev:]
        fnoise_sample = signal[-nstdev:]
        sigerr = fit.noise_stdev(noise_sample)
        fsigerr = fit.noise_stdev(fnoise_sample)
        fitres['noise sigma (mV)'] = sigerr
        amplitude_est = np.mean(signal_unf[:25])
        offset_est = np.mean(signal_unf[-25:])
        frequency_est = fit.estimate_frequency(timebase, signal_unf, 
            name=name)
        fitres['spectral peak'] = frequency_est
        fitres['delta f / f'] = fit.spectrum_fmhw(signal_unf, fs, name=name)
        damping_est = fit.estimate_damping(timebase, signal_unf, 
            frequency_est, name=name)
        # With those estimates ready, we can try fitting the signal
        pguess = [0*amplitude_est, frequency_est, 0*damping_est, 0, 0, 0, 
            0*offset_est, 1]
        lbound=(-100, 0, 0, -100, 0, -pi/2, -100, -100)
        ubound=(100, 6, np.inf, 100, np.inf, pi/2, 100, 100)
        zbounds = list(zip(ubound, lbound))

        if configs.getboolean('fit-params','fit smooth'):
            fitd = signal
        else:
            fitd = signal_unf

        if configs.getboolean('fit-params', 'use shotgun-lsq'):
            spreads = str_to_floats(configs.get('fit-params', 
                'shotgun-lsq spreads'))
            bestp = pguess
            bestp = fit.shotgun_lsq(timebase, signal_unf, sfit,
                                       p0=pguess,
                                       spread=spreads,
                                       sigma=sigerr,
                                       maxiter=1000)
        try:
            bestp, c2r = fit.basin_lsq(sfit, timebase, 
                signal_unf, bestp, sigma=sigerr,
                bounds=zbounds)
            warnings.filterwarnings('ignore')
            bestps, bestcovs = scipy.optimize.curve_fit(sfit, timebase, 
                fitd,
                sigma=sigerr,
                absolute_sigma=True,
                #bounds=(lbound, ubound),
                method='trf',
                #p0=bestp,
                max_nfev=5000)
            c2rs = fit.redchi2(sfit, timebase, signal_unf, bestps, 
                sigma=sigerr)
            r2s = fit.nlcorr(sfit, timebase, signal_unf, bestps)
        except RuntimeError as e:
            print('FAILED S FIT:')
            #print(e)
            c2rs = np.inf
            r2 = -np.inf
            r2s = -np.inf
            bestps = pguess
            bestcovs = np.ones((len(pguess), len(pguess)))
        try:
            bestpe, bestcove = scipy.optimize.curve_fit(expfit, 
                timebase, fitd,
                sigma=sigerr,
                absolute_sigma=True,
                method='trf',
                p0=bestp)
            c2re = fit.redchi2(expfit, timebase, signal_unf, bestpe, 
                sigma=sigerr)
            r2e = fit.nlcorr(expfit, timebase, signal_unf, bestpe)
            #bestcove = np.ones((len(pguess), len(pguess)))
            bestfite = expfit(timebase, *bestpe)
        except RuntimeError as e:
            print('FAILED EXP FIT:')
            #print(e)
            r2 = -np.inf
            r2e = -np.inf
            bestpe = [1, 1, 1, 1, 1, 1, 1]
            bestcove = np.ones((len(pguess), len(pguess)))
        #print('c2rs: {:.4g}, c2re: {:.4g}'.format(c2rs, c2re))
        #if c2re < c2rs:
        if False: # temporarily disabling exp fit
            bestp = [bestpe[0], 0, bestp[2], 0, 
                bestpe[4], bestp[5], bestp[6]]
            bestfit = expfit(timebase, *bestpe)
            bestfitenv = bestfit
            #bestp = [bestpe[0]]
            c2r = c2re
            r2 = r2e
            chi2upper = chi2uppere
            bestcov = bestcove
            fitxt = 'exp fit'
        else:
            bestp = bestps
            bestfit = sfit(timebase, *bestp)
            penv = copy.deepcopy(bestp)
            penv[1] = 0
            bestfitenv = sfit(timebase, *penv)
            chi2upper = chi2uppers
            bestcov = bestcovs
            c2r = c2rs
            r2 = r2s
            fitxt = 'sin fit'

        warnings.resetwarnings()
        #print('BEFORE', bestp)
        #worstint = fit.conf_chi2(sfit, timebase, signal_unf, bestp, 2, sigma=sigerr)
        #print('AFTER', bestp)
        bestpc = fit.conf1d(bestp, bestcov, 2)
        #bestint = [st[1] for st in bestpc]
        #print('DELTA CONFIDENCE INTERVALS: ', np.subtract(bestint, worstint))

        #bestfit = sfit(timebase, *bestp)
        #c2r = fit.redchi2(sfit, timebase, signal_unf, bestp, sigma=sigerr)
        #fit.error_analysis(timebase, signal_unf, bestfit, name=name)
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
        fitres['frequency estimage'] = frequency_est
        fitres['lambda estimate'] = damping_est
        fitres['interference lambda'] = bestpc[4][0]
        fitres['lambda'] = bestpc[2][0]
        fitres['lambda interval'] = bestpc[2][1]
        fitres['amplitude'] = bestpc[0][0]
        fitres['amplitude interval'] = bestpc[0][1]
        fitres['interference amplitude'] = bestpc[3][0]
        fitres['interference damping'] = bestpc[4][0]
        fitres['time delay'] = bestpc[5][0]
        fitres['DC offset'] = bestpc[6][0]
        fitres['chi square'] = c2r
        fitres['best fit SNR'] = fit.snr(signal_unf, bestfit)

        #justexp = fitres['amplitude'] * np.exp(-0.5*fitres['lambda']*timebase)
        #justcos = np.cos(2*pi*fitres['frequency']*timebase + fitres['time delay'])
        # If this signal matches the goodness-of-fit criteria,
        # add it to the fit results dictionary
        if use_for_fit:
            r[name] = fitres
        fitres['used in fit'] = use_for_fit
        rall[name] = fitres

        print('    {:20s}'.format(name), end=' '*2)
        for pk in ['amplitude', 'frequency', 'lambda', 'chi square']:
            print('{:<9.4g}'.format(fitres[pk]), end=' '*2)
        print('')

        tb1 = np.linspace(-1, 5, 1000)
        bestfit1 = sfit(tb1, *bestp)
        plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(timebase, signal_unf, 'b.', label='data')
        plt.plot(timebase, signal, 'r--', label='smoothed data')
        plt.plot(timebase, bestfit, 'g', label=fitxt)

        #plt.plot(timebase, justexp, 'r', label='just exp')
        #plt.plot(timebase, justcos, 'k', label='just cos')
        #plt.plot(timebase, env(bestfit), 'g--', label=fitxt + ' hilb')
        #plt.plot(timebase, bestfite, 'g', label='exp fit')
        paramstr1 = r'r$^2$ = ' + "{0:.3f}\n".format(r2) + fitxt + '\n'
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
        plt.close(fig)
        del fig, ax
        #pb.update(1)

    # All of the fits are in r, so we'll add that to the analysis object
    analysis.set_fits(r)
    analysis.set_rawfits(rall)
    return None


def sfit(x, a, fp, ldamp, aint, intdamp, t0, v0, a0):
    x = np.array(x)
    y = (a * np.cos(2*pi*fp*x-t0) * np.exp(-0.5*ldamp*x)
         #+ aint*x + v0)
         + aint*np.exp(-0.5*intdamp*x) + a0*x + v0)
    return y

def expfit(x, a, fp, ldamp, aint, intdamp, t0, v0, a0):
    x = np.array(x)
    y =   (a * np.exp(-0.5*ldamp*x)
           + aint*x + v0)
    return y


def odfit(x, a, ldamp, t0, v0):
    x = np.array(x)
    y = a * np.exp(-0.5*abs(ldamp)*(x-t0)) + v0


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

def env(x):
    return np.abs(scipy.signal.hilbert(x))
