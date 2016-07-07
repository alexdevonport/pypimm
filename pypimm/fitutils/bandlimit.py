__author__ = 'alex'
import numpy as np
from math import pi, sqrt


def bandlimit(signal, fs, bw):
    """
    Band-limits the given signal using a filter of
    the given order and frequency cutoffs.

    The filter is a linear-phase Gaussian FIR filter
    with a Hamming window.

    :param signal: signal to be filtered, numpy array-like
    :param fs: sampling frequency, float
    :param bwl: pass band lower limit, float
    :param blu: pass band upper limit, float
    :return: filtered signal, numpy array-like
    """

    # the bandpass filter can be constructed from simpler filters:
    #     Bandpass = lowpass (high cutoff) - lowpass(low cutoff)

    # generate the two cutoff lowpass filters
    omega = 2 * pi * bw / fs
    #flt = ideal_lpf(order, omega)
    flt = gaussian(fs, bw)
    l = len(flt)
    # Apply Hamming window
    #hamwin = hamming(order)
    #flt = np.multiply(flt, hamwin)

    # Convolve our signal with the filter
    filtered = np.convolve(signal, flt)

    # Since the filter is linear phase, it has experienced and
    # order/2-point delay. But other than that, there's no
    # distortion!
    ho = int(l/2)
    filtered = filtered[ho:-(l-ho)+1]

    return filtered


def ideal_lpf(l, corner=1):
    """
    The "ideal" FIR low-pass filter is
    y(n) = wc / pi * sinc(wc * n),
    where wc is the cutoff frequency. This function
    returns the first l terms of this function centered
    about n=0 (i.e. from -l/2 to l/2).
    :param l:
    :param corner:
    :return:
    """
    lh = int(l / 2)
    x = np.array(range(-lh, lh + 1))
    r = corner / pi * np.sinc(x * corner / pi)
    if l % 2 == 0:
        r = np.delete(r, l / 2)
    return r

def gaussian(fs, fc):
    """
    For a given sampling frequency and cutoff frequency,
    returns a tuple containing a gaussian filter six
    standard deviations in length.
    :param fs:
    :param fc:
    :return:
    """
    sig = fs / (2 * pi * fc)
    l = 6 * sig #
    lh = int(l / 2)
    x = np.array(range(-lh, lh + 1))
    r = 1 / (sqrt(2 * pi) * sig) * np.exp(-(x * x) / (2 * sig * sig))
    if l % 2 == 0:
        r = np.delete(r, l / 2)
    return r


def hamming(n):
    """
    Returns a length-n Hamming window.
    :param n:
    :return:
    """
    x = np.array(range(n))
    hm = 0.54 - 0.46 * np.cos(2 * pi * x / (n - 1))
    return hm
