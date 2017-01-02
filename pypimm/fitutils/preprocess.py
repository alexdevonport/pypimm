import math
import statsmodels.api as sm
import scipy.signal
import numpy as np
import fitutils as fit
from math import pi

from fitutils import bandlimit

__author__ = 'alex'

def preprocess(timebase, signal, configs, name=''):
    """

    :param timebase:
    :param signal:
    :return:
    """
    timebase = np.array(timebase)

    # convert the signal to mV
    # Originally, we'd normalize the signal to 1, but I decided
    # I'd rather let the fit take care of amplitude, as I might want
    # that information down the line.
    signal = signal * 1E3

    # convert the timebase to ns
    # originally, use of s and ns for the timbase was inconsistent in this program.
    # It's ns all the way after this.
    timebase = timebase * 1E9
    tunf = timebase
    ts = (timebase[1] - timebase[0])  # sampling frequency
    fs = 1 / ts

    # pre-process data for fitting.
    # The data must be filtered to
    # remove noise. Since we aren't real-time, we may as well band-limit
    # via FFT

    bw = configs.getfloat('preprocess', 'bandlimit')
    unf = signal
    raw = signal
    #signal = bandlimit(signal, fs, bw)
    try:
        #signal = fit.optimal_filter(signal, fs, nullsig, [10, 1, 1], noise, [1, 1, 1], name=name)
        signal = scipy.signal.savgol_filter(signal, 15, 4)
    except RuntimeError:
        signal = bandlimit(signal, fs, bw)

    # skip some data. set data[0] to be the first maximum after the skip, so that the
    # signal of interest resembles a damped cosine. Chose damped cosine over damped sine
    # for ease of normalization.
    # In this stage, we also center the signal about zero (remove the mean, that is).

    # arbitrary skip
    # get values from the config file
    #configfp = os.path.join(configdir, 'pypimmconfig.txt')
    tskip = configs.getfloat('preprocess', 'initial data skip')
    tskip1 = configs.getfloat('preprocess', 'skip after max')
    ttrunc = configs.getfloat('preprocess', 'data fit length')
    tzero = configs.getfloat('preprocess', 'zero mean length')
    nzero = math.ceil(tzero / ts)

    nskip = math.ceil(tskip / ts)
    #nskip = 0
    #nskip1 = math.ceil(tskip1 / ts)
    nskip1 = 0
    ntrunc = min(math.ceil(ttrunc / ts), len(signal))  # we don't want the chosen length to be longer than the signal!
    signal = signal[nskip:]
    unf = unf[nskip:]
    timebase = timebase[nskip:] - timebase[nskip]
    tunf = tunf - tunf[nskip]

    if configs.getboolean('preprocess', 'use global max'):
        pk = globalmax(np.abs(signal))
    else:
        pk = localmax(np.abs(signal))
    tpeak = pk *ts
    #totalskip = pk
    #totalskip = pk + nskip1
    totalskip = 0 # No skipping. MADNESSS
    signal = signal[totalskip:]
    unf = unf[totalskip:]
    timebase = timebase[totalskip:] - timebase[totalskip]



    # truncate the signal
    timebase = timebase[:ntrunc]
    signal = signal[:ntrunc]
    unf = unf[:ntrunc]

    # DC shift the signal so that the end data is centered about zero
    dcs = np.mean(signal[-nzero:])
    signal -= dcs
    dcsu = np.mean(unf[-nzero:])
    unf -= dcsu

    return (timebase, signal, tunf, unf, raw, tpeak)




def localmax(y):
    """
     Finds the index of the first local maximum in the signal y.
    :param y:
    :return:
    """
    for k in range(1,len(y)-1):
        if y[k] > y[k-1] and y[k] > y[k+1]:
            return k
    return None

def globalmax(y):
    bestmax = 0
    bestmaxk = 0
    for k in range(1,len(y)-1):
        if y[k] > y[k-1] and y[k] > y[k+1]:
            if y[k] > bestmax:
                bestmax = y[k]
                bestmaxk = k
    return bestmaxk

def noise(f, p0, p1, p2):
    return p0*np.exp(-f*p1) + p2

def lorentz(x, a, x0, gamma):
    n = a * 0.5 * gamma
    d = np.power(x - x0, 2) + (0.5 * gamma)**2
    return 1 / pi * (n / d)

def nullsig(x, a, x0, gamma):
    return np.zeros(np.size(x))

def lorentz2(x, a0, x00, gamma0, a1, x01, gamma1):
    n0 = a0 * 0.5 * gamma0
    d0 = np.power(x - x00, 2) + (0.5 * gamma0)**2
    n1 = a1 * 0.5 * gamma1
    d1 = np.power(x - x01, 2) + (0.5 * gamma1)**2
    return 1 / pi * (n0 / d0) + 1 / pi * (n1 / d1)

def gaussian(x, a, x0, sigma):
    return a / np.sqrt(2*3.1415) * np.exp(-(x - x0)**2 / (2*sigma))