__author__ = 'alex'
from scipy.optimize import curve_fit
from math import pi
import scipy.stats
import scipy.signal
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
    tfit = 0.75
    nfit = int(tfit*fs)
    env = abs(hilbert(signal))
    #lenv = np.log(env)
    #lsig = np.log(np.abs(signal))
    lsin = np.log(np.abs(np.cos(2*pi*fguess*timebase)))
    #lsig = np.log(np.abs(signal))
    signal = signal - np.mean(signal[:-10])
    lsig = (envelope(signal))
    #env=signal
    #env = bandlimit(env, fs, 2.5E9)

    # normalize the envelope and set the last value to zero, so we can focus
    # just on the time constant.

    # The form of the envelope should be close to y=exp(-g*t), where
    # g is the damping we're after. We can now approximate the
    # damping using a least-squares fit.
    try:
        #slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(timebase[:nfit],lsig[:nfit])
        popt, pcov = curve_fit(expdamp, timebase[30:int(len(timebase)/2)], signal[30:int(len(timebase)/2)])
        damping_estimate1 = popt[0]
        damping_estimate2 = popt[1]
        bestfit = expdamp(timebase, *popt)
        #bestfit = expdamp(timebase, popt[0], popt[1], popt[2])
        #damping_estimate = -2 * slope
    except:
        bestfit = expdamp(timebase, 4, 4, 1, 1, 0)
        damping_estimate1 = 4
        damping_estimate2 = 4

    env = env / env[0]
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(timebase, lsig, timebase, bestfit)
    #plt.plot(timebase, lsin)
    paramstr1 = 'lambda estimate 1: {:.4g}\n'.format(damping_estimate1)
    paramstr2 = 'lambda estimate 2: {:.4g}'.format(damping_estimate2)
    plt.text(0.75, 0.25, paramstr1 + paramstr2,
             horizontalalignment='center',
             verticalalignment='center',
             transform=ax.transAxes)
    plt.savefig(r'./envelopes/'+name+'_env.png')
    plt.clf()
    plt.close(fig)
    del fig, ax
    return 0.5*(damping_estimate1 + damping_estimate2)


def hilbert(x):
    npts = len(x)
    hnpts = int(npts / 2)
    xfft = np.fft.fft(x)
    xfft[0:hnpts] *= 2
    xfft[hnpts:] *= 0
    xh = np.fft.ifft(xfft)
    return xh

def envelope(x):
    return np.abs(scipy.signal.hilbert(x))

def expdamp(x, ldamp1, ldamp2, a1, a2, v0):
    return (a1 * np.exp(-0.5*ldamp1 * x)
          + a2 * np.exp(-0.5*ldamp2 * x) + v0)

def pseudoenvelope(x, window=3):
    xr = np.zeros(np.size(x))
    for k, xk in enumerate(x):
        try:
            xr[k] = np.max(x[k-window:k+window])
        except IndexError as e:
            print(e)
            xr[k] = xk
    return xr