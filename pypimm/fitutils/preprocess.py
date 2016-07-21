import math
import statsmodels.api as sm
import scipy.signal
import numpy as np

from fitutils import bandlimit

__author__ = 'alex'

def preprocess(timebase, signal, configs):
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
    signal = bandlimit(signal, fs, bw)
    #signal = scipy.signal.savgol_filter(signal, 21, 3)
    #signal = sm.nonparametric.lowess(signal, timebase, frac=0.02)[:,1]

    # skip some data. set data[0] to be the first maximum after the skip, so that the
    # signal of interest resembles a damped cosine. Chose damped cosine over damped sine
    # for ease of normalization.
    # In this stage, we also center the signal about zero (remove the mean, that is).

    # arbitrary skip
    # get values from the config file
    #configfp = os.path.join(configdir, 'pypimmconfig.txt')
    tskip = configs.getfloat('preprocess', 'initial data skip')
    ttrunc = configs.getfloat('preprocess', 'data fit length')
    tzero = configs.getfloat('preprocess', 'zero mean length')
    nzero = math.ceil(tzero / ts)

    nskip = math.ceil(tskip / ts)
    ntrunc = min(math.ceil(ttrunc / ts), len(signal))  # we don't want the chosen length to be longer than the signal!
    signal = signal[nskip:]
    unf = unf[nskip:]
    timebase = timebase[nskip:] - timebase[nskip]
    tunf = tunf - tunf[nskip]
    # next zero crossing
    #if configs.getboolean('preprocess', 'use global max'):
    #    pk = globalmax(signal)
    #else:
    #    pk = localmax(signal)
    #signal = signal[pk:]
    #timebase = timebase[pk:] - timebase[pk]

    # truncate the signal
    timebase = timebase[:ntrunc]
    signal = signal[:ntrunc]
    unf = unf[:ntrunc]

    # DC shift the signal so that the end data is centered about zero
    dcs = np.mean(signal[-nzero:])
    signal -= dcs
    unf -= dcs

    return (timebase, signal, tunf, unf, raw)




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