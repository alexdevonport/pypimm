__author__ = 'alex'

from bandlimit import bandlimit
import matplotlib.pyplot as plt
import numpy as np
from math import  pi, sqrt
import time

def risetime(signal, fs=None):
    """
    Returns the 10-90 rise time of signal with
    respect to a given sampling frequency (or the
    number of samples if no frequency is given).

    """
    signal = signal/np.max(signal)
    t1 = np.argmin(abs(signal-0.1))
    t2 = np.argmin(abs(signal-0.9))
    tr = 1 / fs * (t2 - t1)
    return tr

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

def digitalh(bcoeffs, acoeffs):
    omega = np.linspace(0, pi, 1000)
    num = np.zeros(1000) + 0j
    den = np.ones(1000) + 0j
    for k, b in enumerate(bcoeffs):
        num += b * np.exp(1j*k*omega)
    for k, a in enumerate(acoeffs):
        den -= a * np.exp(1j*k*omega)
    r = np.divide(num, den)
    return omega, r


t = np.linspace(0, 10, 1000)
ts = t[1] - t[0]
fs = 1 / ts
x = np.concatenate([np.zeros(100), np.ones(900)])

y = bandlimit(x, fs, 0.5)

print(len(x))
print(len(y))
plt.plot(t, x, label='before')
plt.plot(t, y, label='after')
plt.xlabel('time (ns)')
plt.ylabel('signal (V)')
plt.legend()
plt.savefig(r'./fitlertest.png')
plt.clf()

w, h = digitalh(gaussian(1, 0.4), [])

plt.plot(w, abs(h))
plt.savefig(r'./bodetest.png')
plt.clf()

for k in range(10):
    print('\rOn iteration'+str(k)+'...',end='')
    time.sleep(1)