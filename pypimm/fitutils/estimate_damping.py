__author__ = 'alex'
from scipy.optimize import curve_fit
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt


def estimate_damping(timebase, signal, fguess, name=None):
    """

    :param timebase:
    :param signal:
    :return:
    """

    # TODO: band-limit the signal to a low frequency (< 1GHZ)

    # calculate envelope of signal
    # The envelope of a signal can be calculated as the magnitude
    # of the analytic signal of which our original signal is the real
    # part. We can construct the analytic signal using a discrete Hilbert
    # transform.
    fs = 1 / (timebase[1] - timebase[0])
    nfit = int(fs)
    env = abs(hilbert(signal))
    lenv = np.log(env)
    #env=signal
    #env = bandlimit(env, fs, 2.5E9)

    # normalize the envelope and set the last value to zero, so we can focus
    # just on the time constant.

    # The form of the envelope should be close to y=exp(-g*t), where
    # g is the damping we're after. We can now approximate the
    # damping using a least-squares fit.
    try:
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(timebase[:nfit],lenv[:nfit])
        #popt, pcov = curve_fit(expdamp, timebase[30:int(len(timebase)/2)], env[30:int(len(timebase)/2)], p0=[-0.2, 1.0, 0.0])
        bestfit = timebase * slope + intercept
        #bestfit = expdamp(timebase, popt[0], popt[1], popt[2])
        damping_estimate = slope
    except:
        bestfit = expdamp(timebase, -0.2, 1.0, 0.0)
        damping_estimate = 0.2

    env = env / env[0]
    plt.clf()
    plt.plot(timebase, lenv, timebase, bestfit)
    plt.savefig(r'./envelopes/'+name+'_env.png')
    return damping_estimate


def hilbert(x):
    npts = len(x)
    hnpts = int(npts / 2)
    xfft = np.fft.fft(x)
    xfft[0:hnpts] *= 2
    xfft[hnpts:] *= 0
    xh = np.fft.ifft(xfft)
    return xh


def expdamp(x, p0, p1, p2):
    return p1 * np.exp(p0 * x) + p2
