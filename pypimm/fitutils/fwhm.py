__author__ = 'alex'
import numpy
import scipy.signal
import os
import matplotlib.pyplot as plt

def spectrum_fmhw(signal, fs, name=None):
    """
    Finds the full-width half-max of the largtest peak
    in the power spectrum of the input signal.
    :param signal: Numpy array-like of floats
    :param fs: float
    :return: float
    """
    flist, pspect = scipy.signal.periodogram(signal, fs)
    hm = 0.5 * numpy.max(pspect)
    pspectmh = numpy.subtract(pspect,  hm)
    hmpts = zeros(flist, pspectmh)
    try:
        hml, hmr = hmpts[-2],  hmpts[-1]
        hmw = hmr - hml
    except IndexError:
        hmw = numpy.inf
    maxf = numpy.argmax(pspect)
    dff = hmw / maxf  # delta f over f

    return dff

def zeros(xs, ys):
    zs = []
    yprev = ys[0]
    xprev = xs[0]
    for x, y in zip(xs, ys):
        if y * yprev < 0:
            zs.append(0.5 * (x + xprev))
        yprev = y
        xprev = x
    return zs

def main():
    t = numpy.linspace(0, 10, 10000)
    sn = numpy.sin(t)
    zlist = zeros(t, sn)
    print(zlist)
    return 0

if __name__ == '__main__':
    main()