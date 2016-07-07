import math
import sys
import os
import numpy as np
from get_config_value import get_config_value
from bandlimit import bandlimit
__author__ = 'alex'

def pimm_preprocess_data(timebase, signal, configs):
    """

    :param timebase:
    :param signal:
    :return:
    """

    # pre-process data for fitting.
    # The data must be filtered to
    # remove noise. Since we aren't real-time, we may as well band-limit
    # via FFT

    ts = (timebase[1] - timebase[0])  # sampling frequency
    fs = 1 / ts
    bw = 1E9 * configs['bandlimit']
    signal = bandlimit(signal, fs, bw)

    # skip some data. set data[0] to be the first maximum after the skip, so that the
    # signal of interest resembles a damped cosine. Chose damped cosine over damped sine
    # for ease of normalization.
    # In this stage, we also center the signal about zero (remove the mean, that is).
    timebase = np.array(timebase)
    # arbitrary skip
    # get values from the config file
    #configfp = os.path.join(configdir, 'pypimmconfig.txt')
    tskip = 1E-9 * configs['initial data skip']
    ttrunc = 1E-9 * configs['data fit length']

    nskip = math.ceil(tskip / ts)
    ntrunc = min(math.ceil(ttrunc / ts), len(signal))  # we don't want the chosen length to be longer than the signal!
    timebase = timebase[nskip:] - timebase[nskip]
    skipmean = np.mean(signal[nskip:])
    signal = signal[nskip:] - skipmean
    # next zero crossing
    if configs['use global max']:
        pk = globalmax(signal)
    else:
        pk = localmax(signal)
    peakmean = np.mean(signal[pk:])
    timebase = timebase[pk:] - timebase[pk]
    signal = signal[pk:] - peakmean
    # truncate the signal
    timebase = timebase[:ntrunc]
    signal = signal[:ntrunc]
    signal = signal - np.mean(signal)
    # convert the signal to mV
    # Originally, we'd normalize the signal to 1, but I decided
    # I'd rather let the fit take care of amplitude, as I might want
    # that information down the line.
    signal = signal * 1E3

    # convert the timebase to ns
    # originally, use of s and ns for the timbase was inconsistent in this program.
    # It's ns all the way after this.
    timebase = timebase * 1E9

    return (timebase, signal)

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