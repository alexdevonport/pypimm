import numpy as np
from math import log2, ceil
import matplotlib.pyplot as plt
import scipy.signal
import os
__author__ = 'alex'


def pimm_calculate_frequency(timebase, signal, resolution=0.01, name=None):
    """

    :param resolution:
    :param name:
    :param timebase:
    :param signal:
    :return:
    """

    # zero-pad the array.
    # The resolution of an N-point FFT is Fs/N. We want
    # the resolution to be AT LEAST as good as the specified
    # value, so we take N as the smallest power of 2 that
    # satisfies Fs/N < resolution.

    # since we expect the timebase points to be equally spaced, the difference
    # of any two successive values will be the sampling period
    #resolution *= 1E9  # convert resolution to GHz
    fs = 1 / (timebase[1] - timebase[0])
    npts = 2 ** ceil(log2(fs / resolution))

    # now, we make the signal be npts long. While we're at it, we'll put the values into
    if (len(signal) < npts):
        diff = npts - len(signal)
        signal = np.concatenate((np.array(signal), np.zeros(diff)))
    else:  # for a low enough resolution, the original signal length could be longer than npts
        signal = np.array(signal[0:npts])

    # FFT the padded array and get the power spectrum
    #sfft = np.fft.fft(signal)
    #pspect = abs(sfft)
    #nyq = int(len(pspect)/2)
    #flist = np.array(range(len(pspect))) * fs / npts
    flist, pspect = scipy.signal.periodogram(signal, fs)
    plt.clf()
    #plt.plot(flist[:nyq] , pspect[:nyq])
    plt.plot(flist, pspect)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('PSD (V$^2$ / GHz)')
    plt.title(name + ' Spectrum')
    fp = os.path.join('.','spectra', name + 'spectrum.png')
    plt.savefig(fp)
    # find the maximum spectrum value & its index
    # since we're looking at real data, we'll only consider the first half of the spectrum
    maxind = (pspect[:int(len(pspect)/2)]).argmax()
    maxf = fs / npts * maxind
    # return frequency bin corresponding to the max index
    return maxf
